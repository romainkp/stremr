## TODO: Create a table of counterfactual predictions (indexed by t,d,V)

# ---------------------------------------------------------------------------------------
#' Estimate any parametric MSM with previously fitted weights.
#'
#' Estimate the causal MSM for a particular stochastic, dynamic or static intervention on the
#' treatment/exposure and monitoring processes based on
#' the user-specified formula. 
#'
#' @param wts_data A list of \code{data.table}s, each data set is a result of calling the function
#' \code{getIPWeights}. Must contain the treatment/monitoring rule-specific estimated IPTW weights.
#' This argument can be also a single \code{data.table} obtained with \code{data.table::rbindlist(wts_data)}.
#' @param form Formula for the logistic MSM
#' @param family MSM family link function, passed directly to glm().
#' @param OData The object returned by function \code{fitPropensity}. Contains the input dat and the
#' previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param use_weights Logical value. Set to \code{FALSE} to ignore the weights in \code{wts_data} and
#' fit a "crude" MSM that does not adjust for the possible confounding due to non-random
#' assignment of the exposure/censoring and monitoring indicators.
#' @param stabilize Set to \code{TRUE} to stabilize the weights by the empirical conditional probability
#' of having followed the rule at time-point \code{t}, given the subject has followed the rule all the way up to
#' time-point \code{t}.
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the
#' value in \code{trunc_weights} will be truncated.
#' @param weights Optional \code{data.table} with additional observation-time-specific weights.
#' Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param getSEs A logical indicator. Set to \code{TRUE} to evaluate the standard errors for the
#' estimated survival by using the MSM influence curve.
# @param glm_package Which R package should be used for fitting the weighted logistic regression
# model (MSM) for the survival hazard?
# Currently available options are \code{"glm"} only.
#' @param return_wts Return the data.table with subject-specific IP weights as part of the output.
#' Note: for large datasets setting this to \code{TRUE} may lead to extremely large object sizes!
#' @param tvals Vector of time-points for Y(t), default (NULL) is to use all available time-points.
#' @param tmax Maximum value of the follow-up period.
#' All person-time observations above this value will be excluded from the MSM model.
#' @param rule.var The column name in wts_data that identifies each unique regimen, can be
#' integer, continuous, character or factor.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#'
#' @return MSM estimation results
#' @seealso \code{\link{fitPropensity}}, \code{\link{getIPWeights}}.
#' @example tests/examples/4_survMSM_example.R
#' @export
fit_gMSM <- function(wts_data,
                     form,
                     family,
                     OData,
                     use_weights = TRUE,
                     stabilize = TRUE,
                     trunc_weights = 10^6,
                     weights = NULL,
                     getSEs = TRUE,
                     return_wts = FALSE,
                     tvals = NULL,
                     tmax = NULL,
                     rule.var = "rule.name",
                     verbose = getOption("stremr.verbose")) {

  if (missing(form)) stop("must specify MSM 'formula'")
  form <- formula(form)
  if (!(inherits(form, "formula"))) {
    stop("formula(form) must return a formula")
  }

  if (missing(wts_data)) stop("must specify weights dataset as 'wts_data'")
  # glm_package <- glm_package[1L]
  # if (!(glm_package %in% c("glm", "speedglm", "h2o"))) stop("glm_package must be either 'speedglm' or 'h2o'")

  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  ## rbind all regimen-specific weights datasets into single dataset (will be used to form design matrix):
  wts_data <- format_wts_data(wts_data)

  ## keep all weights, even if they are 0 (for correct inference, glm will not use 0 weights)
  ## remove all observations with missing weights or outcomes
  MSM_fit_data <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))), ]
  setkeyv(MSM_fit_data, cols = c(nodes$IDnode, nodes$tnode))

  # multiply the weight by stabilization factor (numerator) (doesn't do anything for saturated MSMs, since cum.stab.P cancels out):
  if (stabilize) MSM_fit_data[, "cum.IPAW" := cum.stab.P * cum.IPAW]

  # truncation of the weights:
  if (trunc_weights < Inf) MSM_fit_data[cum.IPAW > trunc_weights, cum.IPAW := trunc_weights]

  # Add additional (user-supplied) observation-specific weights to the cumulative weights:
  MSM_fit_data <- process_opt_wts(MSM_fit_data, weights, nodes, adjust_outcome = FALSE)

  # When !use_weights run a crude estimator by setting all non-zero weights to 1 and removing all zero weights:
  if (!use_weights) MSM_fit_data[cum.IPAW > 0, cum.IPAW := 1L]

  # Evaluate all time-points observed in the data (t's). Subset by tmax or select specified tvals only
  if (!is.null(tmax)) {
    MSM_fit_data <- MSM_fit_data[eval(as.name(nodes$tnode)) <= tmax, ]
  }

  if (!is.null(tvals)) {
    MSM_fit_data <- MSM_fit_data[eval(as.name(nodes$tnode)) %in% tvals, ]
  }

  t_periods_used <- sort(unique(MSM_fit_data[[t_name]]))

  ## subset weights data only by the time-points being considered:
  if (verbose) { cat("using MSM tvals: "); cat(t_periods_used); cat("\n") }

  ## MSM formula without the outcome ("Y ~ ...")  
  form_woutY <- getCovariateFormula(form)

  ## Vector of all vars used in MSM formula (including V):
  MSMform_vars <- all.vars(as.formula(form_woutY))

  ## Additional vars used in MSM formula that are not t or d (e.g., summary f(t,d))
  new_MSM_vars <- setdiff(MSMform_vars, c(rule.var, t_name))
  if (verbose) message(paste0("new covariates (not t or d): ", paste0(new_MSM_vars, collapse=",")))

  ## Additional vars used in MSM formula that are not part of the weights
  MSM_V_vars <- setdiff(MSMform_vars, names(MSM_fit_data))
  if (verbose) message(paste0("bsl covariates (V): ", paste0(MSM_V_vars, collapse=",")))


  ## TODO 1: Select only on V (MSM_V_vars) that already exist in OData:
  ##         Some Vs are false-positives (e.g, value-named of existing vars as in I(rule.name==gTI.dhigh)))
  ## --------Decided not to do this, most of the time above indicates error in formula spec's

  ## Select V values @ min(t) in OData and join on ID with MSM_fit_data
  tmin <- min(wts_data[[nodes$tnode]], na.rm = TRUE)
  V_data <- OData$dat.sVar[eval(as.name(nodes$tnode)) == tmin, .SD, .SDcols = c(nodes$IDnode, MSM_V_vars)]
  ## Left join MSM_fit_data & V_data with (keeps all rows in MSM_fit_data, even if V value is NA):
  MSM_fit_data <- V_data[MSM_fit_data,, on = nodes$IDnode]
  setcolorder(MSM_fit_data, c(setdiff(names(MSM_fit_data), MSM_V_vars), MSM_V_vars))

  ## direct MSM approach with arbitrary link fun
  ctrl <- stats::glm.control(trace = FALSE)
  msm.fit <- glm(form, family = family, data = MSM_fit_data, control = ctrl, weights = cum.IPAW)
  msm.fit$data <- NULL
  design_mat <- model.matrix(as.formula(form_woutY), data = MSM_fit_data)

  # if (verbose) message("...evaluating MSM-based survival curves...")

  ## obtain hazard predictions for each observation O_i:
  MSM_fit_data[, "glm.IPAW.predict" := predict(msm.fit, newdata = MSM_fit_data, type = "response")]

  ## unique values of the rule / regimen identifier:
  # rules_TRT <- sort(unique(MSM_fit_data[[rule.var]]))
  # dt_names <- list(t_periods_used, rules_TRT)
  # names(dt_names) <- c(t_name, rule.var)
  # estimates <- data.table(purrr::cross_df(dt_names))
  # setkeyv(estimates, c(rule.var, t_name))

  ## todo: Need to add a check that there is no more than one unique value of f(d,t) for each d and t
  # MSM_vars_by_t_d <- wts_data[, unique(.SD), by = c(rule.var, t_name), .SDcols = new_MSM_vars]
  # ## add additional MSM covariates/summaries to the by-rule prediction matrix:
  # if (length(new_MSM_vars) > 0) {
  #   estimates <- estimates[MSM_vars_by_t_d, ]
  # }

  # ## obtain hazard predictions for all (d,t) combos:
  # estimates[, "ht.MSM" := predict(msm.fit, newdata = estimates, type = "response")]
  # ## evaluate survival from hz(d,t) for each regimen in d:
  # estimates[, "St.MSM" := cumprod(1-ht.MSM), by = eval(as.name(rule.var))]

  #### For variance (SEs), GET IC and SE FOR BETA's
  #### GET IC and SE FOR Sd(t)
  ## S.d.t.predict - MSM survival estimates for one regimen
  ## h.d.t.predict - MSM hazard estimates for one regimen
  ## design.d.t - d-specific matrix of dummy indicators for each t, i.e., d(m(t,d))/t
  ## IC.O - observation-specific IC estimates for MSM coefs
  ## IC.S - observation-specific IC estimates for S(t) (by time-point)

  ## this is a matrix of partial derivatives of the regression formula (wrt each beta_j)
  ## needed for IC of survival (St)
  ## each row is indexing some time-point (t) and regimen (d)
  # new_design_mat <- model.matrix(as.formula(form_woutY), data = estimates)

  if (getSEs) {
    if (verbose) message("...evaluating SEs based on MSM coefficients fits and the estimated IC...")
    # IC.Var.S.d <- vector(mode = "list", length(rules_TRT))
    # names(IC.Var.S.d) <- rules_TRT

    ## the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
    ## the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
    ## Always use IPW_MSMestimator=TRUE (even for crude estimators) to only use rule-followers.
    beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                               MSMdata = MSM_fit_data,
                               MSMdesign = design_mat,
                               MSMpredict = "glm.IPAW.predict",
                               IPW_MSMestimator = TRUE)

    # for(d.j in rules_TRT) {
    #   rule.idx <- estimates[[rule.var]] %in% d.j
    #   IC.Var.S.d[[d.j]] <- getSE.S(nID = nID,
    #                                S.d.t.predict = estimates[rule.idx, St.MSM],
    #                                h.d.t.predict = estimates[rule.idx, ht.MSM],
    #                                design.d.t = new_design_mat[rule.idx,],
    #                                IC.O = beta.IC.O.SEs[["IC.O"]])

    #   estimates[estimates[[rule.var]]==d.j, ("SE.MSM") := IC.Var.S.d[[d.j]][["se.S"]]]
    #   n_ts <- nrow(IC.Var.S.d[[d.j]][["IC.S"]])
    #   est_dj <- estimates[estimates[[rule.var]]==d.j, ]

    #   for (i in 1:n_ts) {
    #     est_dj[i, ("IC.St") := list(list(IC.Var.S.d[[d.j]][["IC.S"]][i, ]))]
    #   }
    #   estimates[estimates[[rule.var]]==d.j, ("IC.St") := est_dj[["IC.St"]]]
    # }
  } else {
    beta.IC.O.SEs <- NULL
  }

  ## The output MSM object.
  ##  "estimates" is a data.table with surv estimates (column "St.MSM");
  ##  Separate row for each time-point t;
  ##  "estimates" contains IC.St (observation-specific IC estimates for S(t)) saved in a column
  # est_name <- "MSM"
  # attr(estimates, "estimator_short") <- est_name
  # attr(estimates, "estimator_long") <- "MSM for E(Y_{d,t}|V)"
  # attr(estimates, "nID") <- nID
  # # attr(estimates, "time") <- estimates[[t_name]]
  # attr(estimates, "trunc_weights") <- trunc_weights

  # return(list(estimates = estimates, msm.fit = msm.fit, beta.SE = beta.IC.O.SEs$se.beta))
  return(list(msm.fit = msm.fit, beta.SE = beta.IC.O.SEs$se.beta))

}
