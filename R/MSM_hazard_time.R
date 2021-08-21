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
#' Estimate Survival with a particular MSM for the survival-hazard function using previously fitted weights.
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the
#' treatment/exposure and monitoring processes based on
#' the user-specified Marginal Structural Model (MSM) for the counterfactual survival function.
#'
#' @param wts_data A list of \code{data.table}s, each data set is a result of calling the function
#' \code{getIPWeights}. Must contain the treatment/monitoring rule-specific estimated IPTW weights.
#' This argument can be also a single \code{data.table} obtained with \code{data.table::rbindlist(wts_data)}.
#' @param OData The object returned by function \code{fitPropensity}. Contains the input dat and the
#' previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param tbreaks The vector of integer (or numeric) breaks that defines the dummy indicators of the
#' follow-up in the observed data. Used for fitting the parametric (or saturated) MSM for
#' the survival hazard. See Details.
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
#' Currently available options are \code{"glm"}, \code{"speedglm"} and \code{"h2o"}.
#' \code{h2o} can provided better performance
#' when fitting MSM with many observations and large number of time-points.
#' @param return_wts Return the data.table with subject-specific IP weights as part of the output.
#' Note: for large datasets setting this to \code{TRUE} may lead to extremely large object sizes!
#' @param tmax Maximum value of the follow-up period.
#' All person-time observations above this value will be excluded from the MSM model.
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
survMSM <- function(wts_data,
                    OData,
                    tbreaks,
                    use_weights = TRUE,
                    stabilize = TRUE,
                    trunc_weights = 10^6,
                    weights = NULL,
                    getSEs = TRUE,
                    est_name = "IPW",
                    glm_package = c("glm", "speedglm", "h2o"),
                    return_wts = FALSE,
                    tmax = NULL,
                    verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  wts_data <- format_wts_data(wts_data)
  rules_TRT <- sort(unique(wts_data[["rule.name"]]))

  glm_package <- glm_package[1L]
  if (!(glm_package %in% c("glm", "speedglm", "h2o"))) stop("glm_package must be either 'glm', 'speedglm' or 'h2o'")

  if (verbose) print("performing MSM estimation for the following TRT/MONITOR rules found in column 'rule.name': " %+% paste(rules_TRT, collapse=","))

  ## Remove all observations with 0 cumulative weights & copy the weights data.table
  ## keep all weights, even if they are 0:
  wts_data_used <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))), ]
  ## remove cumulative weights that are 0:
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

  # Default tbreaks, error checks for tbreaks, plus padding w/ mint & tmax:
  if (missing(tbreaks)) {
    # default tbreaks is to use a saturated (non-parametric) MSM
    tbreaks <- sort(periods)
    if (verbose)
      message("running MSM with default 'tbreaks': (" %+%
        paste0(tbreaks, collapse = ",") %+%
        "); \nNote: such 'tbreaks' define a separate coefficient for every unique follow-up time period resulting in a saturated (non-parametric) MSM.")
  }

  if (length(unique(tbreaks)) < length(tbreaks))
    stop("all tbreaks must be unique")

  if (!all(tbreaks %in% periods))
    stop("all tbreaks must be contained between minimum and maximum follow-up periods:" %+% tbreaks[!(tbreaks %in% periods)])

  if (max(tbreaks) < tmax) tbreaks <- sort(c(tbreaks, tmax)) # pad on the right (if needed with tmax):

  # Create the dummies I(d == intervened_TRT) for the logistic MSM for d-specific hazard
  all.d.dummies <- NULL
  for( dummy.j in rules_TRT ){
    wts_data_used[, (dummy.j) := as.integer(rule.name %in% dummy.j)]
    all.d.dummies <- c(all.d.dummies, dummy.j)
  }

  # Create the dummies I(t in interval.j), where interval.j defined by intervals of time of increasing length
  all.t.dummies <- NULL
  tbreaks.mint <- c(mint, tbreaks) # pad tbreaks on the left (with mint)
  MSM.intervals <- matrix(NA, ncol = 2, nrow = length(tbreaks)) # save the actual intervals
  colnames(MSM.intervals) <- c("min.t", "max.t")
  t.per.inteval <- vector(mode = "list", length = nrow(MSM.intervals)) # save the vector of period vals that belong to each interval
  for (t.idx in 2:length(tbreaks.mint)) {
    low.t <- tbreaks.mint[t.idx - 1]
    high.t <- tbreaks.mint[t.idx]
    # First interval needs to be closed on both sides (includes the smallest obesrved follow-up, mint)
    if (t.idx == 2L) {
      dummy.j <- paste("Periods.", low.t, "to", high.t, sep="")
      MSM.intervals[t.idx - 1, ] <- c(low.t, high.t); t.per.inteval[[t.idx - 1]] <- unique(low.t:high.t)
      wts_data_used[, (dummy.j) := as.integer(get(t_name) >= low.t & get(t_name) <= high.t)]
      # wts_data_used[, (dummy.j) := as.integer(eval(as.name(t_name)) >= low.t & eval(as.name(t_name)) <= high.t)]
    } else {
      dummy.j <- paste("Periods.", (low.t + 1), "to", high.t, sep="")
      MSM.intervals[t.idx - 1, ] <- c(low.t + 1, high.t)
      t.per.inteval[[t.idx - 1]] <- unique((low.t+1):high.t)
      wts_data_used[, (dummy.j) := as.integer(get(t_name) >= (low.t + 1) & get(t_name) <= high.t)]
      # wts_data_used[, (dummy.j) := as.integer(eval(as.name(t_name)) >= (low.t + 1) & eval(as.name(t_name)) <= high.t)]
    }
    if (verbose) print("defined t.dummy: " %+% dummy.j)
    all.t.dummies <- c(all.t.dummies, dummy.j)
  }

  # Create interaction dummies I(t in interval.j & d == intervened_TRT)
  for (d.dummy in all.d.dummies) {
    for (t.dummy in all.t.dummies) {
      if (verbose) print(t.dummy %+% "_" %+% d.dummy)
      wts_data_used[, (t.dummy %+% "_" %+% d.dummy) := as.integer(eval(as.name(t.dummy)) & eval(as.name(d.dummy)))]
    }
  }

  all_dummies <-  paste(sapply(all.d.dummies, function(x) {
                        return(paste(paste(paste(all.t.dummies, x, sep="_"), sep="")))
                        }))

  # Fit the hazard MSM
  resglmMSM <- runglmMSM(wts_data_used, all_dummies, Ynode, glm_package, verbose)
  wts_data_used[, glm.IPAW.predictP1 := resglmMSM$glm.IPAW.predictP1]
  m.fit <- resglmMSM$m.fit

  # Compute the Survival curves under each d
  S2.IPAW <- hazard.IPAW <- rep(list(rep(NA,tmax-mint+1)), length(rules_TRT))
  names(S2.IPAW) <- names(hazard.IPAW) <- rules_TRT

  if (verbose) message("...evaluating MSM-based survival curves...")

  for(d.j in names(S2.IPAW)) {
    for(period.idx in seq_along(periods)) {
      period.j <- periods[period.idx] # the period of the follow-up for which we want to evaluate the MSM-based survival:
      # the dummy coefficient of the MSM that includes this time-point (period)
      # that is, find the appropriate right-closed interval from MSM.intervals matrix for a given period.j:
      int_idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
      if (!(period.j %in% t.per.inteval[[int_idx]])) stop("error while finding the appropriate MSM coefficient")
      d.j.idx <- which(all.d.dummies %in% d.j)
      MSM.term <- all_dummies[length(all.t.dummies)*(d.j.idx - 1) + int_idx]
      # print("fetching MSM coefficient for period " %+% period.j %+% " and rule " %+% d.j %+% ": " %+% MSM.term)
      hazard.IPAW[[d.j]][period.idx] <- 1 / (1 + exp(-m.fit$coef[MSM.term]))
      S2.IPAW[[d.j]][period.idx] <- (1-1/(1 + exp(-m.fit$coef[MSM.term])))
    }
    S2.IPAW[[d.j]] <- cumprod(S2.IPAW[[d.j]])
  }

  if (verbose) message("...evaluating SEs based on MSM hazard fit and the estimated IC...")

  #### For variance (SEs), GET IC and SE FOR BETA's
  #### GET IC and SE FOR Sd(t)
  # S.d.t.predict - MSM survival estimates for one regimen
  # h.d.t.predict - MSM hazard estimates for one regimen
  # design.d.t - d-specific matrix of dummy indicators for each t, i.e., d(m(t,d))/t
  # IC.O - observation-specific IC estimates for MSM coefs
  # IC.S - observation-specific IC estimates for S(t) (by time-point)
  design.d.t <- rep(list(matrix(0L, ncol = length(all_dummies), nrow = length(mint:tmax))), length(rules_TRT))
  IC.Var.S.d <- vector(mode = "list", length(rules_TRT))
  names(design.d.t) <- names(IC.Var.S.d) <- rules_TRT

  for(d.j in names(S2.IPAW)) {
    IC.Var.S.d[[d.j]][["se.S"]] <- rep(NA_real_, length(S2.IPAW[[d.j]]))
    IC.Var.S.d[[d.j]][["IC.S"]] <- matrix(NA_real_, nrow = length(periods), ncol = nID)
  }

  ## wrapping the whole thing in try({}), if it fails, use empty se/IC above
  if (getSEs) {
    try({
      # the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
      # the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
      # p.coef - number of time-specific coefficients in the MSM
      p.coef <- nrow(MSM.intervals) # p.coef <- length(tjmin)
      design.t <- matrix(0L, ncol = p.coef, nrow = length(periods))
      for (period.idx in seq_along(periods)) {
        period.j <- periods[period.idx]
        col.idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
        design.t[period.idx, col.idx] <- 1
      }
      
      ## Always use IPW_MSMestimator=TRUE (even for crude estimators) to only use rule-followers.
      beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                                MSMdata = wts_data_used, MSMdesign = as.matrix(wts_data_used[, all_dummies, with = FALSE]),
                                MSMpredict = "glm.IPAW.predictP1", IPW_MSMestimator = TRUE)

      for(d.j in names(S2.IPAW)) {
        d.idx <- which(names(S2.IPAW) %in% d.j)
        set_cols <- seq((d.idx - 1) * ncol(design.t) + 1, (d.idx) * ncol(design.t))
        design.d.t[[d.j]][,set_cols] <- design.t

        IC.Var.S.d[[d.j]] <- getSE.S(nID = nID,
                                     S.d.t.predict = S2.IPAW[[d.j]],
                                     h.d.t.predict = hazard.IPAW[[d.j]],
                                     design.d.t = design.d.t[[d.j]],
                                     IC.O = beta.IC.O.SEs[["IC.O"]])
      }      
    })
  }

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
      for (i in 1:n_ts) {
        estimates[i, ("IC.St") := list(list(IC.Var.S.d[[rule_name]][["IC.S"]][i, ]))]
      }

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
  if (length(MSM_out) == 1L) return(MSM_out[[1L]]) else  return(MSM_out)
}

runglmMSM <- function(wts_data, all_dummies, Ynode, glm_package, verbose) {
  # Generic prediction fun for logistic regression coefs, predicts P(A = 1 | X_mat)
  # Does not handle cases with deterministic Anodes in the original data.
  logispredict = function(m.fit, X_mat) {
    eta <- X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
    pAout <- match.fun(FUN = m.fit$linkfun)(eta)
    return(pAout)
  }

  if (glm_package %in% "h2o") {
    if (verbose) message("...fitting hazard MSM with h2o::h2o.glm...")
    MSM.designmat.H2O <- fast.load.to.H2O(wts_data, destination_frame = "MSM.designmat.H2O")
    m.fit_h2o <- try(h2o::h2o.glm(y = Ynode,
                                  x = all_dummies,
                                  intercept = FALSE,
                                  weights_column = "cum.IPAW",
                                  training_frame = MSM.designmat.H2O,
                                  family = "binomial",
                                  standardize = FALSE,
                                  solver = "IRLSM", # solver = c("L_BFGS"),
                                  lambda = 0L,
                                  max_iterations = 50,
                                  ignore_const_cols = FALSE
                                  ),
                silent = TRUE)

    out_coef <- vector(mode = "numeric", length = length(all_dummies))
    out_coef[] <- NA
    names(out_coef) <- c(all_dummies)
    out_coef[names(m.fit_h2o@model$coefficients)[-1]] <- m.fit_h2o@model$coefficients[-1]
    m.fit <- list(coef = out_coef, linkfun = "plogis", fitfunname = "h2o.glm")
    glm.IPAW.predictP1 <- as.vector(h2o::h2o.predict(m.fit_h2o, newdata = MSM.designmat.H2O)[,"p1"])
    # wts_data[, glm.IPAW.predictP1 := as.vector(h2o::h2o.predict(m.fit_h2o, newdata = MSM.designmat.H2O)[,"p1"])]

  } else if (glm_package %in% c("speedglm","glm")) {

    if (verbose) message(paste0("...fitting hazard MSM with ", glm_package))
    Xdesign.mat <- as.matrix(wts_data[, all_dummies, with = FALSE])

    if (glm_package %in% "speedglm") {
      m.fit <- try(speedglm::speedglm.wfit(
                                         X = Xdesign.mat,
                                         y = as.numeric(wts_data[[Ynode]]),
                                         intercept = FALSE,
                                         family = quasibinomial(),
                                         weights = wts_data[["cum.IPAW"]],
                                         trace = FALSE),
                          silent = TRUE)
      if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
        if (verbose) stop("speedglm::speedglm.wfit failed, try using glm_package='glm'", m.fit)
        return(invisible(m.fit))
      }      
    }

    if (glm_package %in% "glm") {
      ctrl <- stats::glm.control(trace = FALSE)
      SuppressGivenWarnings({
        m.fit <- stats::glm.fit(x = Xdesign.mat,
                                y = as.numeric(wts_data[[Ynode]]),
                                family = quasibinomial(),
                                intercept = FALSE,
                                weights = wts_data[["cum.IPAW"]],
                                control = ctrl)
      }, GetWarningsToSuppress())
    }

    m.fit <- list(coef = m.fit$coef, linkfun = "plogis", fitfunname = "speedglm")
    if (verbose) {print("MSM fits"); print(m.fit$coef)}
    glm.IPAW.predictP1 <- logispredict(m.fit, Xdesign.mat)
    # wts_data[, glm.IPAW.predictP1 := logispredict(m.fit, Xdesign.mat)]
  }

  return(list(glm.IPAW.predictP1 = glm.IPAW.predictP1, m.fit = m.fit))
}