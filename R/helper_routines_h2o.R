# nocov start

# helper function for h2o frames
h2o.plogis <- function(x) {
  h2o_exp_x <- h2o::h2o.exp(x)
  h2o_exp_x / (1 + h2o_exp_x)
}

# helper function for h2o frames
h2o.qlogis <- function(x) h2o::h2o.log(x / (1 - x))

# Poor mans slicer until h2o error is fixed:
reassign_rows_cols <- function(data, newsubset_idx, subset_frame, col_name = "Qkplus1", ID = "ID", t = "t") {
  # full rows in h2o frame that will need to be set:
  copy_rows_to_set <- data$H2Oframe[newsubset_idx, ]
  # remove those rows from main frame:
  data$H2Oframe <- data$H2Oframe[-newsubset_idx,]
  # set the colum in subset frame:
  copy_rows_to_set[, "Qkplus1"] <- subset_frame
  # add the subset frame back into the frame with rbind
  data$H2Oframe <- h2o::h2o.rbind(data$H2Oframe, copy_rows_to_set)
  # sort by t:
  data$H2Oframe <- h2o::h2o.arrange(data$H2Oframe, ID, t)
  return(invisible(NULL))
}

#************************************************
# TMLE update with H2OFrame and h2o.glm
#************************************************
tmle.update.h2o <- function(prev_Qkplus1, init_Q_fitted_only, IPWts, lower_bound_zero_Q = TRUE, skip_update_zero_Q = TRUE) {
  QY.star <- NA
  # h2o.sum(h2o.abs(h2oFRAME_test[, "test_col_copy"]))
  if (sum(abs(IPWts)) < 10^-9) {
  # if (h2o.sum(h2o.abs(IPWts)) < 10^-9) {
    update.Qstar.coef <- 0
    if (gvars$verbose) message("TMLE update cannot be performed since all IP-weights are exactly zero!")
    warning("TMLE update cannot be performed since all IP-weights are exactly zero!")
  } else if ((sum(prev_Qkplus1[IPWts > 0]) < 10^-5) && skip_update_zero_Q) {
  # } else if ((h2o.sum(prev_Qkplus1[IPWts > 0]) < 10^-5) && skip_update_zero_Q) {
    update.Qstar.coef <- 0
  } else {
    #************************************************
    # TMLE update via weighted univariate ML (espsilon is intercept)
    #************************************************
    if (lower_bound_zero_Q) {
      prev_Qkplus1[prev_Qkplus1 < 10^(-4)] <- 10^(-4)
      init_Q_fitted_only[init_Q_fitted_only < 10^(-4)] <- 10^(-4)
    }

    off_TMLE <- h2o.qlogis(init_Q_fitted_only)

    # Create a column (h2o frame) with consant vector of 1
    # fit_h2oframe <- h2o::h2o.rep_len(1L, nrow(prev_Qkplus1))
    # fit_h2oframe <- h2o::h2o.cbind(fit_h2oframe, prev_Qkplus1, off_TMLE, IPWts)
    # names(fit_h2oframe) <- c("Intercept", "y", "offset_column", "weights_column")

    fit_h2oframe <- h2o::h2o.cbind(prev_Qkplus1, off_TMLE, IPWts)
    # fit_h2oframe <- cbind(prev_Qkplus1, off_TMLE, IPWts)
    names(fit_h2oframe) <- c("y", "offset_column", "weights_column")
    # x = "Intercept",
    # browser()

    # solver = c("AUTO",
    #   "IRLSM", "L_BFGS", "COORDINATE_DESCENT_NAIVE", "COORDINATE_DESCENT")
    LFM_fit_t <- system.time(
    m.Qstar <- h2o::h2o.glm(y = "y", training_frame = fit_h2oframe,
                            lambda = 0L, ignore_const_cols = FALSE, intercept = TRUE, standardize = FALSE,
                            family = "quasibinomial", offset_column = "offset_column", weights_column = "weights_column")
      )
    print("LFM_fit_t"); print(LFM_fit_t)

    if (inherits(m.Qstar, "try-error")) { # TMLE update failed
      if (gvars$verbose) message("attempt at running TMLE update with speedglm::speedglm.wfit has failed")
      warning("attempt at running TMLE update with speedglm::speedglm.wfit has failed")
      update.Qstar.coef <- 0
    } else {
      update.Qstar.coef <- m.Qstar@model$coefficients
    }
  }

  fit <- list(TMLE.intercept = update.Qstar.coef)
  class(fit)[2] <- "tmlefit"
  # if (gvars$verbose)
    print("tmle update: " %+% update.Qstar.coef)
  return(fit)
}

# nocov end