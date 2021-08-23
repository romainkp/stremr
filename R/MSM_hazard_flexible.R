if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("ht.MSM", "St.MSM"))
}

#' @importFrom stats formula gaussian model.matrix
NULL

getCovariateFormula <- function(object) {
  form <- formula(object)
  if (!(inherits(form, "formula"))) {
    stop("formula(object) must return a formula")
  }
  form <- form[[length(form)]]
  if (length(form) == 3 && form[[1]] == as.name("|")) {
    form <- form[[2]]
  }
  eval(substitute(~form))
}

format_wts_data <- function(wts_data) {
  # Stack the weighted data sets, if those came in a list:
  if (is.data.table(wts_data)) {
    # ...do nothing...
  } else if (is.list(wts_data)) {
    assert_that(all(sapply(wts_data, is.data.table)))
    wts_data <- rbindlist(wts_data)
  } else {
    stop("...wts_data arg must be either a list of rule-specific weight data.tables or one combined weight data.table...")
  }
return(wts_data)
}

# ---------------------------------------------------------------------------------------
#' Estimate Survival with any MSM for the survival-hazard function using previously fitted weights.
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the
#' treatment/exposure and monitoring processes based on
#' the user-specified Marginal Structural Model (MSM) for the counterfactual survival function.
#'
#' @param wts_data A list of \code{data.table}s, each data set is a result of calling the function
#' \code{getIPWeights}. Must contain the treatment/monitoring rule-specific estimated IPTW weights.
#' This argument can be also a single \code{data.table} obtained with \code{data.table::rbindlist(wts_data)}.
#' @param form Formula for the logistic MSM
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
#' @param glm_package Which R package should be used for fitting the weighted logistic regression
#' model (MSM) for the survival hazard?
#' Currently available options are \code{"glm"} only.
#' @param return_wts Return the data.table with subject-specific IP weights as part of the output.
#' Note: for large datasets setting this to \code{TRUE} may lead to extremely large object sizes!
#' @param tmax Maximum value of the follow-up period.
#' All person-time observations above this value will be excluded from the MSM model.
#' @param rule.var The column name in wts_data that identifies each unique regimen, can be
#' integer, continuous, character or factor.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#'
#' @section MSM:
#' **********************************************************************
#' This routine will run the weighted logistic regression using the (possibly-weighted) outcomes from
#' many regimens, with dummy indicators for each treatment/monitoring
#' regimen available in \code{wts_data} and each follow-up time interval specified in \code{tbreaks}.
#' When \code{use_weights = TRUE}, the logistic regression for the survival hazard is weighted by the
#' \strong{IPW} (Inverse Probability-Weighted or Horvitz-Thompson) estimated weights
#' in \code{wts_data}. These IPW weights are based on the previously fitted propensity scores (function
#' \code{fitPropensity}), allowing
#' adjustment for confounding by possibly non-random assignment to exposure and monitoring and possibly
#' informative right-censoring.
#'
#' @section Specifying time-intervals:
#' **********************************************************************
#'
#' \code{tbreaks} is used for defining the time-intervals of the MSM coefficients for estimation of the
#' survival hazard function.
#' The first value in \code{tbreaks} defines a dummy variable (indicator) for a fully closed interval,
#' with each subsequent value in \code{tbreaks} defining a single right-closed time-interval.
#' For example, \code{tbreaks = c(0,1)} will define the MSM dummy indicators: I(\code{tmin} <= t <= 0 ),
#' I(0 < t <= 1) and I(1 < t <= \code{tmax}),
#' where \code{tmin} is the minimum value of the time variable (automatically determined from input weights)
#' and \code{tmax} is the maximum value of the time variable ( if omitted this will also be automatically
#' determined from the input weights).
#'
#' On the other hand \code{tbreaks = c(1)} will define the following (more parametric) MSM dummy
#' indicators: I(\code{mint} <= t <=1 ) and I(1 < t <= \code{tmax}).
#' If omitted, the default is to define a saturated (non-parametric) MSM with a separate dummy variable
#' for every unique period in the observed data.
#'
#' @return MSM estimation results composed of a separate list for each treatment regimen.
#' Each regimen-specific list contains an item named \code{"estimates"}, which is a data.table
#' with MSM survival estimates in a column \code{"St.MSM"}. The data.table \code{"estimates"} contains
#' a separate row for each time-point \code{t}. The \code{"estimates"} also contains the
#' standard error estimates for MSM survival and the observation-specific influence-curve estimates for
#' the MSM survival saved in a column named \code{"IC.St"}.
#' @seealso \code{\link{fitPropensity}}, \code{\link{getIPWeights}}.
#' @example tests/examples/4_survMSM_example.R
#' @export
fit_hMSM <- function(wts_data,
                     form,
                     OData,
                     use_weights = TRUE,
                     stabilize = TRUE,
                     trunc_weights = 10^6,
                     weights = NULL,
                     getSEs = TRUE,
                     glm_package = c("glm", "speedglm", "h2o"),
                     return_wts = FALSE,
                     tmax = NULL,
                     rule.var = "rule.name",
                     verbose = getOption("stremr.verbose")) {

  if (missing(form)) stop("must specify MSM 'formula' for the hazard")
  form <- formula(form)
  if (!(inherits(form, "formula"))) {
    stop("formula(form) must return a formula")
  }
  form_woutY <- getCovariateFormula(form)
  MSMform_vars <- all.vars(as.formula(form_woutY))

  if (missing(wts_data)) stop("must specify weights dataset as 'wts_data'")
  glm_package <- glm_package[1L]
  if (!(glm_package %in% c("glm", "speedglm", "h2o"))) stop("glm_package must be either 'speedglm' or 'h2o'")

  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  ## rbind all regimen-specific weights datasets into one:
  wts_data <- format_wts_data(wts_data)

  ## remove all observations with 0 cumulative weights & copy the weights data.table
  ## keep all weights, even if they are 0:
  wts_data_used <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))), ]
  ## remove cumulative weights that are 0 (can result in an error for saturated MSM):
  # wts_data_used <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))) & (cum.IPAW > 0), ]
  setkeyv(wts_data_used, cols = c(nodes$IDnode, nodes$tnode))

  # multiply the weight by stabilization factor (numerator) (doesn't do anything for saturated MSMs, since cum.stab.P cancels out):
  if (stabilize) wts_data_used[, "cum.IPAW" := cum.stab.P * cum.IPAW]

  # If trunc_weights < Inf, do truncation of the weights
  if (trunc_weights < Inf) wts_data_used[cum.IPAW > trunc_weights, cum.IPAW := trunc_weights]

  # Add additional (user-supplied) observation-specific weights to the cumulative weights:
  wts_data_used <- process_opt_wts(wts_data_used, weights, nodes, adjust_outcome = FALSE)

  # When !use_weights run a crude estimator by setting all non-zero weights to 1
  if (!use_weights) wts_data_used[cum.IPAW > 0, cum.IPAW := 1L]

  # Define all observed sequence of periods (t's)
  mint <- min(wts_data_used[[t_name]], na.rm = TRUE)
  if (is.null(tmax)) tmax <- max(wts_data_used[[nodes$tnode]], na.rm = TRUE)

  ## subset weights data only by the time-points being considered:
  wts_data_used <- wts_data_used[eval(as.name(nodes$tnode)) <= tmax, ]

  ## fit the IPW-MSM for the hazard
  ## direct GLM MSM approach
  ctrl <- stats::glm.control(trace = FALSE)
  msm.fit <- glm(form, family = quasibinomial(), data = wts_data_used, control = ctrl, weights = cum.IPAW)
  # msm.fit <- speedglm::speedglm(form, family = quasibinomial(), data = wts_data_used, weights = wts_data_used[["cum.IPAW"]])
  msm.fit$data <- NULL
  design_mat <- model.matrix(as.formula(form_woutY), data = wts_data_used)

  if (verbose) message("...evaluating MSM-based survival curves...")

  ## obtain hazard predictions for each observation O_i:
  wts_data_used[, "glm.IPAW.predictP1" := predict(msm.fit, newdata = wts_data_used, type = "response")]

  ## create the output data.table with estimates (by time / regimen)
  periods <- (mint:tmax)

  ## unique values of the rule / regimen identifier:
  rules_TRT <- sort(unique(wts_data[[rule.var]]))

  ## evaluate additional vars used in MSM formula, as f(d,t)
  MSMform_vars <- all.vars(as.formula(form_woutY))
  new_MSM_vars <- setdiff(MSMform_vars, c(rule.var, t_name))
  ## todo: Need to add a check that there is no more than one unique value of f(d,t) for each d and t
  MSM_vars_by_t_d <- wts_data[, unique(.SD), by = c(rule.var, t_name), .SDcols = new_MSM_vars]

  if (verbose) { print("MSM t values: "); print(periods) }

  dt_names <- list(periods, rules_TRT)
  names(dt_names) <- c(t_name, rule.var)
  estimates <- data.table(purrr::cross_df(dt_names))
  setkeyv(estimates, c(rule.var, t_name))

  ## add additional MSM covariates/summaries to the by-rule prediction matrix:
  if (length(new_MSM_vars) > 0) {
    estimates <- estimates[MSM_vars_by_t_d, ]
  }

  ## obtain hazard predictions for all (d,t) combos:
  estimates[, "ht.MSM" := predict(msm.fit, newdata = estimates, type = "response")]
  ## evaluate survival from hz(d,t) for each regimen in d:
  estimates[, "St.MSM" := cumprod(1-ht.MSM), by = eval(as.name(rule.var))]

  if (verbose) message("...evaluating SEs based on MSM hazard fit and the estimated IC...")

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
  new_design_mat <- model.matrix(as.formula(form_woutY), data = estimates)

  if (getSEs) {
    IC.Var.S.d <- vector(mode = "list", length(rules_TRT))
    names(IC.Var.S.d) <- rules_TRT

    ## the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
    ## the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
    ## Always use IPW_MSMestimator=TRUE (even for crude estimators) to only use rule-followers.
    beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                              MSMdata = wts_data_used,
                              # MSMdesign = as.matrix(wts_data_used[, all_dummies, with = FALSE]),
                              MSMdesign = design_mat,
                              MSMpredict = "glm.IPAW.predictP1",
                              IPW_MSMestimator = TRUE)

    for(d.j in rules_TRT) {
      rule.idx <- estimates[[rule.var]] %in% d.j
      IC.Var.S.d[[d.j]] <- getSE.S(nID = nID,
                                   S.d.t.predict = estimates[rule.idx, St.MSM],
                                   h.d.t.predict = estimates[rule.idx, ht.MSM],
                                   design.d.t = new_design_mat[rule.idx,],
                                   IC.O = beta.IC.O.SEs[["IC.O"]])

      estimates[estimates[[rule.var]]==d.j, ("SE.MSM") := IC.Var.S.d[[d.j]][["se.S"]]]
      n_ts <- nrow(IC.Var.S.d[[d.j]][["IC.S"]])
      est_dj <- estimates[estimates[[rule.var]]==d.j, ]

      for (i in 1:n_ts) {
        est_dj[i, ("IC.St") := list(list(IC.Var.S.d[[d.j]][["IC.S"]][i, ]))]
      }
      estimates[estimates[[rule.var]]==d.j, ("IC.St") := est_dj[["IC.St"]]]
    }
  } else {
    IC.Var.S.d <- NULL
  }

  ## The output MSM object.
  ##  "estimates" is a data.table with surv estimates (column "St.MSM");
  ##  Separate row for each time-point t;
  ##  "estimates" contains IC.St (observation-specific IC estimates for S(t)) saved in a column
  est_name <- "MSM"
  attr(estimates, "estimator_short") <- est_name
  attr(estimates, "estimator_long") <- "MSM (Marginal Structural Model) for hazard, mapped into survival"
  attr(estimates, "nID") <- nID
  attr(estimates, "time") <- estimates[[t_name]]
  attr(estimates, "trunc_weights") <- trunc_weights

  return(list(estimates = estimates, msm.fit = msm.fit, beta.SE = beta.IC.O.SEs$se.beta))


}
