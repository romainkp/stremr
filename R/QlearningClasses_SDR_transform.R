# nocov start
SDRtransform <- R6Class(classname = "SDRtransform",
  inherit = ModelGeneric,
  portable = TRUE,
  class = TRUE,
  public = list(
    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      # serial loop over all regressions in PsAsW.models:
      max_Qk_idx <- length(private$PsAsW.models)
      ## need to set ignore_tmax, which should be the max follow-up time-point
      for (k_i in seq_along(private$PsAsW.models)) {
        ## 1. Fit the inital for current k_i (Q_k)
        private$PsAsW.models[[k_i]]$fit(data = data, ...)
        ## 2. Perform SDR tranform that will define the outcome for the next regression.
        ## This is done iteratively, over all previous Q_k fits, starting from the very first initial regression (Q_K).
        data$dat.sVar[, ("EIC_i_t_sum") := 0.0]
        ignore_tmin <- private$PsAsW.models[[k_i]]$t_period
        wts <- eval_weights_k(data = data, ignore_tmin = ignore_tmin, reverse_wt_prod = FALSE, ...)

        # SDR_term_k_cum <- 0
        for (i in (1:k_i)) {
            SDR_term_k <- private$PsAsW.models[[i]]$transform_Q_k(data = data,
                                                    k_i = k_i,
                                                    i = i,
                                                    max_Qk_idx = max_Qk_idx,
                                                    ...)
            # SDR_term_k_cum <- SDR_term_k_cum + SDR_term_k
        }
      }
      invisible(self)
    }
  )
)

# Q1: What happens when we have reached the final time-point?

SDRtransformQModel  <- R6Class(classname = "SDRtransformQModel",
  inherit = ModelQlearn,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    ## SDR tranform learner
    ## Transform the most recent fit of Q[k] (Qk_hat) into DR outcome (scaled by the reverse cumulative product of the weights)
    ## Move that outcome to next regression step (subset_idx-1)

    ## Redefine the cumulative sum for current k:
    ## Q.weighted.sum[k] = wts * (Qkplus1 - Qk_hat) + Q.weighted.sum[k-1]
    ## The transform needs to be of the following form:
    ## Qk_hat + Q.weighted.sum[k]
    ## This redefines the Qkplus1 for the next regression to follow.

    ## Outcome: The same outcomes that were used to fit this initial Q[k']
    ## Weights: the product of g's from t=k to current k'
    transform_Q_k = function(data, k_i, i, max_Qk_idx, ...) {
      ## use only the observations that participated in fitting of the initial Q_{k'} (current time-point is k')
      use_subset_idx <- self$idx_used_to_fit_initQ
      if (gvars$verbose == 2) {
        cat("...running SDR targeting loop...\n")
        cat("Total number of time-points to consider: " %+% max_Qk_idx, "\n")
        cat("Current (S)DR k index = " %+% self$all_Qregs_indx[k_i], "\n")
        cat("Current (S)DR k' index = " %+% self$all_Qregs_indx[i], "\n")
      }

      ## 1. Weights: defined new column of cumulative weights where cumulative product starts at t = Qk_idx (k), rather than t = 0:
      # wts <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      wts <- data$IPwts_by_regimen[self$subset_idx, cum.IPAW]
      if (self$reg$stabilize) wts <- wts * data$IPwts_by_regimen[self$subset_idx, cum.stab.P]

      ## 2. Outcome: **TARGETED** prediction of the previous step k'+1.
      # Qkplus1.protected <- data$dat.sVar[use_subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      Qkplus1.protected <- data$dat.sVar[self$subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      ## 3. Initial model fits:
      # Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]
      Qk_hat <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
      ## The cumulative transform sum (making sure we are using Qkplus1 and not Gamma transforms):
      SDR_term_k <- wts * (Qkplus1.protected - Qk_hat)
      # print("mean(SDR_term_k)"); print(mean(SDR_term_k))

      # wts_fitted <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      # Qkplus1.protected_fitted <- data$dat.sVar[use_subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      # Qk_hat_fitted <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]
      # SDR_term_k_fitted <- wts_fitted * (Qkplus1.protected_fitted - Qk_hat_fitted)
      # print("mean(SDR_term_k_fitted)"); print(mean(SDR_term_k_fitted))

      data$dat.sVar[self$subset_idx, ("EIC_i_t_sum") := EIC_i_t_sum + SDR_term_k]

      ## Signifies we are in a final (S)DR transform loop (final / last regression for t=0 has been reached).
      ## No more fits of the initials are going to be performed, just evaluate the EIC
      if (k_i == max_Qk_idx) {
        # data$dat.sVar[use_subset_idx, ("EIC_i_t") := SDR_term_k_fitted]
        data$dat.sVar[self$subset_idx, ("EIC_i_t") := SDR_term_k]
      }

      ## Reached the last iter of current SDR loop.
      ## This is where we actually define (S)DR transform Q^*[k].
      ## This is our current model update and it will be used as the outcome for the next regression.
      if (i == k_i) {
        Gamma_DR <- Qk_hat + data$dat.sVar[self$subset_idx, EIC_i_t_sum]
        # Gamma_DR <- Qk_hat_fitted + data$dat.sVar[use_subset_idx, EIC_i_t_sum]
      }

      if ((self$Qreg_counter > 1) && (i == k_i)) {
        # if (gvars$verbose) cat("updating the outcomes for next regression with DR tranformed outcomes (only among those who were at risk and following treatment)")
        data$dat.sVar[self$subset_idx - 1, "Qkplus1.protected" := Qkplus1]
        data$dat.sVar[self$subset_idx - 1, "Qkplus1" := Gamma_DR]

      } else if (self$Qreg_counter > 1) {
        data$dat.sVar[self$subset_idx - 1, ("EIC_i_t_sum") := data$dat.sVar[self$subset_idx, EIC_i_t_sum]]

      } else {
        data$dat.sVar[self$subset_idx, ("Qk_hat") := Gamma_DR]
        private$probAeqa <- Gamma_DR
      }

      if (i == k_i) {
        # print("mean(Gamma_DR)"); print(mean(Gamma_DR))
        return(Gamma_DR)
      } else {
        return(SDR_term_k)
      }
    }
  )
)
# nocov end