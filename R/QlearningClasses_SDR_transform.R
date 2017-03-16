SDRtransform <- R6Class(classname = "SDRtransform",
  inherit = GenericModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      # serial loop over all regressions in PsAsW.models:
      max_Qk_idx <- length(private$PsAsW.models)
      ## need to set ignore_tmax, which should be the max follow-up time-point
      # private$PsAsW.models[[k_i]]$eval_weights_k(data = data, ignore_tmax = ignore_tmax, reverse_wt_prod = TRUE, ...)

      # wts[1:100, ]
      # browser()

      for (k_i in seq_along(private$PsAsW.models)) {
        # browser()
        ## 1. Fit the inital for current k_i (Q_k)
        private$PsAsW.models[[k_i]]$fit(data = data, ...)

        ## 2. Perform SDR tranform that will define the outcome for the next regression.
        ## This is done iteratively, over all previous Q_k fits, starting from the very first initial regression (Q_K).
        data$dat.sVar[, ("EIC_i_t_sum") := 0.0]
        ignore_tmin <- private$PsAsW.models[[k_i]]$t_period
        wts <- eval_weights_k(data = data, ignore_tmin = ignore_tmin, reverse_wt_prod = FALSE, ...)
        # wts <- eval_weights_k(data = data, ignore_tmax = ignore_tmax, reverse_wt_prod = FALSE, ...)

        for (i in (1:k_i)) {
          # private$PsAsW.models[[i]]$transform_Q_k(data = data,
          #                                           k_i = k_i,
          #                                           max_Qk_idx = max_Qk_idx,
          #                                           ...)
          private$PsAsW.models[[i]]$transform_Q_k(data = data,
                                                  k_i = k_i,
                                                  i = i,
                                                  max_Qk_idx = max_Qk_idx,
                                                  # QModel_h_k = QModel_h_k,
                                                  # QModel_Qkplus1 = QModel_Qkplus1,
                                                  ...)
        }
      }
      invisible(self)
    }
  )
)

# Q1: What happens when we have reached the final time-point?

SDRtransformQModel  <- R6Class(classname = "SDRtransformQModel",
  inherit = QlearnModel,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    split_preds_Qk_hat = NULL,
    SDR_SL_fit = NULL,

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
      # if (gvars$verbose == 2) {
        cat("...running SDR targeting loop...\n")
        cat("Total number of time-points to consider: " %+% max_Qk_idx, "\n")
        cat("Current (S)DR k index = " %+% self$all_Qregs_indx[k_i], "\n")
        cat("Current (S)DR k' index = " %+% self$all_Qregs_indx[i], "\n")
      # }
      # browser()

      data$IPwts_by_regimen[use_subset_idx,][1:100, ]
      data$IPwts_by_regimen[self$subset_idx, ][1:100, ]

      data$dat.sVar[use_subset_idx,][1:100, ]
      data$dat.sVar[self$subset_idx,][1:100, ]

      ## 1. Weights: defined new column of cumulative weights where cumulative product starts at t = Qk_idx (k), rather than t = 0:
      # wts <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      wts <- data$IPwts_by_regimen[self$subset_idx, "cum.IPAW", with = FALSE][[1]]
      ## 2a. Outcome: **TARGETED** prediction of the previous step k'+1.
      # Qkplus1 <- data$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]
      Qkplus1 <- data$dat.sVar[self$subset_idx, "Qkplus1", with = FALSE][[1]]
      ## 3. Initial model fits:
      # Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]
      Qk_hat <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
      ## 3. Initial model fit that includes extrapolated predictions (newsly censored and non-rule followers):
      # Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]

      ## The EIC estimates and the cumulative transform sum:
      ## REPLACE Qk_hat with Qk_hat_transform (targeted version)
      SDR_term_k <- wts * (Qkplus1 - Qk_hat)
      # data$dat.sVar[use_subset_idx, ("EIC_i_t") := SDR_term_k]
      # data$dat.sVar[use_subset_idx, EIC_i_t_sum]
      # data$dat.sVar[use_subset_idx,]

      # data$dat.sVar[use_subset_idx, ("EIC_i_t_sum") := EIC_i_t_sum + SDR_term_k]
      data$dat.sVar[self$subset_idx, ("EIC_i_t_sum") := EIC_i_t_sum + SDR_term_k]

      # data$dat.sVar[self$subset_idx, ]

      ## Signifies we are in a final (S)DR transform loop (final / last regression for t=0 has been reached).
      ## No more fits of the initials are going to be performed, just evaluate the EIC
      if (k_i == max_Qk_idx) {
        # data$dat.sVar[use_subset_idx, ("EIC_i_t") := SDR_term_k]
        data$dat.sVar[self$subset_idx, ("EIC_i_t") := SDR_term_k]
      }
      ## Signifies we reached a final iteration of the final (S)DR transform loop. No more initial fits are going to be performed.
      if (i == max_Qk_idx) {

      }

      ## Reached the last iter of current SDR loop.
      ## This is where we actually define (S)DR transform Q^*[k].
      ## This is our current model update and it will be used as the outcome for the next regression.
      if (i == k_i) {
        # Qk_hat <- Qk_hat + data$dat.sVar[use_subset_idx, EIC_i_t_sum]
        Qk_hat_transform <- Qk_hat + data$dat.sVar[self$subset_idx, EIC_i_t_sum]
        # Qk_hat <- Qk_hat + SDR_term_k
        # data$dat.sVar[use_subset_idx, "Qk_hat" := Qk_hat]
        data$dat.sVar[self$subset_idx, "Qk_hat" := Qk_hat_transform]
      }

      # data$dat.sVar[use_subset_idx,]
      # data$dat.sVar[self$subset_idx,]
      # data$dat.sVar[use_subset_idx - 1,]
      # data$dat.sVar[self$subset_idx - 1,]
      # browser()

      # Over-write the old predictions with new model updates as Qk_hat[k'] in row [k']:
      # data$dat.sVar[self$subset_idx, "Qk_hat" := Qk_hat_star_all]
      # Qk_hat_star_all-data$dat.sVar[self$subset_idx, Qk_hat]
      # print("MSE of previous Qk_hat vs. upated Qk_hat: " %+% mean((Qk_hat_star_all-Qk_hat_all)^2))

      if ((self$Qreg_counter > 1) && (i == k_i)) {
        # if (gvars$verbose)
        # Qk_hat_transform <- Qk_hat + data$dat.sVar[use_subset_idx, EIC_i_t_sum]
        cat("updating the outcomes for next regression with DR tranformed outcomes (only among those who were at risk and following treatment)")
        # data$dat.sVar[use_subset_idx - 1, "Qkplus1" := Qk_hat]
        # data$dat.sVar[self$subset_idx - 1, "Qkplus1" := data$dat.sVar[self$subset_idx, "Qk_hat"]]
        data$dat.sVar[self$subset_idx - 1, "Qkplus1" := Qk_hat_transform]

      } else if (self$Qreg_counter > 1) {
        ## Q: Is this right????
        ## I.e., is it possible to have subjects that previously contributed to rolling sum, but have no contributed
        ## at this time-point?
        # data$dat.sVar[self$subset_idx - 1, ("Qk_hat") := Qk_hat_transform]
        data$dat.sVar[self$subset_idx - 1, ("EIC_i_t_sum") := data$dat.sVar[self$subset_idx, EIC_i_t_sum]]

      } else {
        # data$dat.sVar[use_subset_idx, ("Qk_hat") := Qk_hat]
        # data$dat.sVar[self$subset_idx, ("Qk_hat") := Qk_hat_transform]
        private$probAeqa <- data$dat.sVar[self$subset_idx, Qk_hat]
      }

      invisible(self)
    }
  )
)