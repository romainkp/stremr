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
#' @param est_name A string naming the current MSM estimator. Ignored by the current routine but is
#' used when generating reports with \code{make_report_rmd}.
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
survMSM2 <- function(wts_data,
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


  if (missing(wts_data)) stop("must specify weights dataset as 'wts_data'")
  glm_package <- glm_package[1L]
  if (!(glm_package %in% c("glm", "speedglm", "h2o"))) stop("glm_package must be either 'speedglm' or 'h2o'")

  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  wts_data <- format_wts_data(wts_data)
  rules_TRT <- sort(unique(wts_data[[rule.var]]))

  ## Remove all observations with 0 cumulative weights & copy the weights data.table
  ## keep all weights, even if they are 0:
  wts_data_used <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))), ]
  ## remove cumulative weights that are 0 (can result in an error for saturated MSM):
  # wts_data_used <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))) & (cum.IPAW > 0), ]
  setkeyv(wts_data_used, cols = c(nodes$IDnode, nodes$tnode))

  # Multiply the weight by stabilization factor (numerator) (doesn't do anything for saturated MSMs, since cum.stab.P cancels out):
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

  periods <- (mint:tmax)
  periods_idx <- seq_along(periods)
  if (verbose) { print("MSM periods: "); print(periods) }

  ## direct GLM MSM approach
  ctrl <- stats::glm.control(trace = FALSE)
  msm.fit <- glm(form, family = quasibinomial(), data = wts_data_used, weights = cum.IPAW, control = ctrl)
  # msm.fit <- speedglm::speedglm(form, family = quasibinomial(), data = wts_data_used, weights = wts_data_used[["cum.IPAW"]])

  msm.fit$data <- NULL
  # str(msm.fit)
  # msm.fit$terms
  # head(msm.fit$model)
  design_mat <- model.matrix(as.formula(form_woutY), data = wts_data_used)
  # head(design_mat)

  ## Fit the hazard MSM
  # resglmMSM <- runglmMSM(wts_data_used, all_dummies, Ynode, glm_package, verbose)
  # wts_data_used[, glm.IPAW.predictP1 := resglmMSM$glm.IPAW.predictP1]
  # m.fit <- resglmMSM$m.fit
  ## Compute the Survival curves under each d
  # S2.IPAW <- hazard.IPAW <- rep(list(rep(NA,tmax-mint+1)), length(rules_TRT))
  # names(S2.IPAW) <- names(hazard.IPAW) <- rules_TRT

  if (verbose) message("...evaluating MSM-based survival curves...")

  # for(d.j in names(S2.IPAW)) {
  #   for(period.idx in seq_along(periods)) {
  #     period.j <- periods[period.idx] # the period of the follow-up for which we want to evaluate the MSM-based survival:
  #     # the dummy coefficient of the MSM that includes this time-point (period)
  #     # that is, find the appropriate right-closed interval from MSM.intervals matrix for a given period.j:
  #     int_idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
  #     if (!(period.j %in% t.per.inteval[[int_idx]])) stop("error while finding the appropriate MSM coefficient")
  #     d.j.idx <- which(all.d.dummies %in% d.j)
  #     MSM.term <- all_dummies[length(all.t.dummies)*(d.j.idx - 1) + int_idx]
  #     # print("fetching MSM coefficient for period " %+% period.j %+% " and rule " %+% d.j %+% ": " %+% MSM.term)
  #     hazard.IPAW[[d.j]][period.idx] <- 1 / (1 + exp(-m.fit$coef[MSM.term]))
  #     S2.IPAW[[d.j]][period.idx] <- (1-1/(1 + exp(-m.fit$coef[MSM.term])))
  #   }
  #   S2.IPAW[[d.j]] <- cumprod(S2.IPAW[[d.j]])
  # }

  ## obtain hazard predictions for each observation O_i:
  wts_data_used[, "glm.IPAW.predictP1" := predict(msm.fit, newdata = wts_data_used, type = "response")]

  dt_names <- list(periods, rules_TRT)
  names(dt_names) <- c(t_name, rule.var)
  new_dat <- data.table(purrr::cross_d(dt_names))
  setkeyv(new_dat, c(rule.var, t_name))

  ## obtain hazard predictions for all (d,t) combos:
  new_dat[, "haz" := predict(msm.fit, newdata = new_dat, type = "response")]
  ## evaluate survival from hazard for each regimen in d:
  new_dat[, "St" := cumprod(1-haz), by = eval(as.name(rule.var))]

  if (verbose) message("...evaluating SEs based on MSM hazard fit and the estimated IC...")

  #### For variance (SEs), GET IC and SE FOR BETA's
  #### GET IC and SE FOR Sd(t)
  # S.d.t.predict - MSM survival estimates for one regimen
  # h.d.t.predict - MSM hazard estimates for one regimen
  # design.d.t - d-specific matrix of dummy indicators for each t, i.e., d(m(t,d))/t
  # IC.O - observation-specific IC estimates for MSM coefs
  # IC.S - observation-specific IC estimates for S(t) (by time-point)


  ## this is a matrix of partial derivatives of the regression formula (wrt each beta_j)
  ## needed for IC of survival (St)
  ## each row is indexing some time-point (t) and regimen (d)
  new_design_mat <- model.matrix(as.formula(form_woutY), data = new_dat)

  if (getSEs) {
    # design.d.t <- rep(list(matrix(0L, ncol = length(all_dummies), nrow = length(mint:tmax))), length(rules_TRT))
    IC.Var.S.d <- vector(mode = "list", length(rules_TRT))
    # names(design.d.t) <- names(IC.Var.S.d) <- rules_TRT

    # # the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
    # # the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
    # # p.coef - number of time-specific coefficients in the MSM
    # p.coef <- nrow(MSM.intervals) # p.coef <- length(tjmin)
    # design.t <- matrix(0L, ncol = p.coef, nrow = length(periods))
    # for (period.idx in seq_along(periods)) {
    #   period.j <- periods[period.idx]
    #   col.idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
    #   design.t[period.idx, col.idx] <- 1
    # }

    beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                              MSMdata = wts_data_used,
                              # MSMdesign = as.matrix(wts_data_used[, all_dummies, with = FALSE]),
                              MSMdesign = design_mat,
                              MSMpredict = "glm.IPAW.predictP1",
                              IPW_MSMestimator = use_weights)

    for(d.j in rules_TRT) {
      rule.idx <- new_dat[[rule.var]] %in% d.j
      # d.idx <- which(names(S2.IPAW) %in% d.j)
      # set_cols <- seq((d.idx - 1) * ncol(design.t) + 1, (d.idx) * ncol(design.t))
      # design.d.t[[d.j]][,set_cols] <- design.t

      IC.Var.S.d[[d.j]] <- getSE.S(nID = nID,
                                   S.d.t.predict = new_dat[rule.idx, St],
                                   h.d.t.predict = new_dat[rule.idx, haz],
                                   design.d.t = new_design_mat[rule.idx,],
                                   IC.O = beta.IC.O.SEs[["IC.O"]])
    }
  } else {
    IC.Var.S.d <- NULL
  }

  return(list(new_dat, IC.Var.S.d))
  browser()

  # # Default tbreaks, error checks for tbreaks, plus padding w/ mint & tmax:
  # if (missing(tbreaks)) {
  #   # default tbreaks is to use a saturated (non-parametric) MSM
  #   tbreaks <- sort(periods)
  #   if (verbose)
  #     message("running MSM with default 'tbreaks': (" %+%
  #       paste0(tbreaks, collapse = ",") %+%
  #       "); \nNote: such 'tbreaks' define a separate coefficient for every unique follow-up time period resulting in a saturated (non-parametric) MSM.")
  # }

  # if (length(unique(tbreaks)) < length(tbreaks))
  #   stop("all tbreaks must be unique")

  # if (!all(tbreaks %in% periods))
  #   stop("all tbreaks must be contained between minimum and maximum follow-up periods:" %+% tbreaks[!(tbreaks %in% periods)])

  # if (max(tbreaks) < tmax) tbreaks <- sort(c(tbreaks, tmax)) # pad on the right (if needed with tmax):

  # # Create the dummies I(d == intervened_TRT) for the logistic MSM for d-specific hazard
  # all.d.dummies <- NULL
  # for( dummy.j in rules_TRT ){
  #   wts_data_used[, (dummy.j) := as.integer(rule.name %in% dummy.j)]
  #   all.d.dummies <- c(all.d.dummies, dummy.j)
  # }

  # # Create the dummies I(t in interval.j), where interval.j defined by intervals of time of increasing length
  # all.t.dummies <- NULL
  # tbreaks.mint <- c(mint, tbreaks) # pad tbreaks on the left (with mint)
  # MSM.intervals <- matrix(NA, ncol = 2, nrow = length(tbreaks)) # save the actual intervals
  # colnames(MSM.intervals) <- c("min.t", "max.t")
  # t.per.inteval <- vector(mode = "list", length = nrow(MSM.intervals)) # save the vector of period vals that belong to each interval
  # for (t.idx in 2:length(tbreaks.mint)) {
  #   low.t <- tbreaks.mint[t.idx - 1]
  #   high.t <- tbreaks.mint[t.idx]
  #   # First interval needs to be closed on both sides (includes the smallest obesrved follow-up, mint)
  #   if (t.idx == 2L) {
  #     dummy.j <- paste("Periods.", low.t, "to", high.t, sep="")
  #     MSM.intervals[t.idx - 1, ] <- c(low.t, high.t); t.per.inteval[[t.idx - 1]] <- unique(low.t:high.t)
  #     wts_data_used[, (dummy.j) := as.integer(get(t_name) >= low.t & get(t_name) <= high.t)]
  #     # wts_data_used[, (dummy.j) := as.integer(eval(as.name(t_name)) >= low.t & eval(as.name(t_name)) <= high.t)]
  #   } else {
  #     dummy.j <- paste("Periods.", (low.t + 1), "to", high.t, sep="")
  #     MSM.intervals[t.idx - 1, ] <- c(low.t + 1, high.t)
  #     t.per.inteval[[t.idx - 1]] <- unique((low.t+1):high.t)
  #     wts_data_used[, (dummy.j) := as.integer(get(t_name) >= (low.t + 1) & get(t_name) <= high.t)]
  #     # wts_data_used[, (dummy.j) := as.integer(eval(as.name(t_name)) >= (low.t + 1) & eval(as.name(t_name)) <= high.t)]
  #   }
  #   if (verbose) print("defined t.dummy: " %+% dummy.j)
  #   all.t.dummies <- c(all.t.dummies, dummy.j)
  # }

  # # Create interaction dummies I(t in interval.j & d == intervened_TRT)
  # for (d.dummy in all.d.dummies) {
  #   for (t.dummy in all.t.dummies) {
  #     if (verbose) print(t.dummy %+% "_" %+% d.dummy)
  #     wts_data_used[, (t.dummy %+% "_" %+% d.dummy) := as.integer(eval(as.name(t.dummy)) & eval(as.name(d.dummy)))]
  #   }
  # }

  # all_dummies <-  paste(sapply(all.d.dummies, function(x) {
  #                       return(paste(paste(paste(all.t.dummies, x, sep="_"), sep="")))
  #                       }))

  ## The output MSM object.
  ## "estimates" is a data.table with surv estimates (column "St.MSM");
  ##  Separate row for each time-point t;
  ##  "estimates" contains IC.St (observation-specific IC estimates for S(t)) saved in a column
  MSM_out <- lapply(rules_TRT, function(rule_name) {
      est_name <- "MSM"
      estimates <- data.table(time = periods,
                              ht.MSM = hazard.IPAW[[rule_name]],
                              St.MSM = S2.IPAW[[rule_name]],
                              SE.MSM = IC.Var.S.d[[rule_name]][["se.S"]],
                              rule.name = rep(rule_name, length(periods)))
      estimates <- cbind(est_name = est_name, estimates)

      n_ts <- nrow(IC.Var.S.d[[rule_name]][["IC.S"]])
      for (i in 1:n_ts)
        estimates[i, ("IC.St") := list(list(IC.Var.S.d[[rule_name]][["IC.S"]][i, ]))]

      attr(estimates, "estimator_short") <- est_name
      attr(estimates, "estimator_long") <- "MSM (Marginal Structural Model) for hazard, mapped into survival"
      attr(estimates, "nID") <- nID
      attr(estimates, "rule_name") <- rule_name
      attr(estimates, "time") <- estimates[["time"]]
      attr(estimates, "trunc_weights") <- trunc_weights

      return(list(est_name = est_name,
                  periods = periods,
                  St = S2.IPAW[[rule_name]],
                  ht = hazard.IPAW[[rule_name]],
                  MSM.fit = m.fit,
                  MSM.intervals = MSM.intervals,
                  # IC.Var.S.d = IC.Var.S.d[[rule_name]],
                  nID = nID,
                  nobs = nrow(wts_data_used),
                  wts_data = { if (return_wts) {wts_data_used} else {NULL} },
                  use_weights = use_weights,
                  trunc_weights = trunc_weights,
                  estimates = estimates
                  ))
    }
  )

  names(MSM_out) <- rules_TRT

  if (length(MSM_out) == 1L) return(MSM_out[[1L]]) else return(MSM_out)
}
