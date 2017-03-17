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

        # if (k_i < max_Qk_idx) {
        SDR_term_k_cum <- 0
        for (i in (1:k_i)) {
            # Qk_hat <- data$dat.sVar[t==1, Qk_hat]
            # Qk_hat <- data$dat.sVar[(private$PsAsW.models[[2]]$idx_used_to_fit_initQ), Qk_hat]
            # Qk_hat <- data$dat.sVar[(private$PsAsW.models[[2]]$idx_used_to_fit_initQ-1), Qkplus1]
            # Qk_hat <- data$dat.sVar[(private$PsAsW.models[[1]]$idx_used_to_fit_initQ), Qk_hat]
            # mean(Qk_hat)
            # wts <- wts[private$PsAsW.models[[2]]$idx_used_to_fit_initQ, cum.IPAW]
            # Qkplus1 <- data$dat.sVar[private$PsAsW.models[[2]]$idx_used_to_fit_initQ, Qkplus1]
            # Qk_hat <- data$dat.sVar[t==0, Qk_hat]
            # mean(Qk_hat)
            # mean(data$dat.sVar[t==1, Qkplus1])
            # mean(wts*(Qkplus1-Qk_hat))

            # private$PsAsW.models[[i]]$idx_used_to_fit_initQ
            # (wts[t==1, cum.IPAW] * (data$dat.sVar[t==1, Qkplus1] - Qk_hat))
            # mean(wts[t==1, cum.IPAW] * (data$dat.sVar[t==1, Qkplus1] - Qk_hat))
            # # mean(SDR_term_k)

            # mean((wts[t==1, cum.IPAW] * (data$dat.sVar[t==1, Qkplus1] - Qk_hat))[wts[t==1, cum.IPAW] > 0])
            # mean

            # data$dat.sVar[t==1, Qk_hat := data$dat.sVar[t==0, Qkplus1]]
            # mean(Qk_hat + SDR_term_k)
            # mean(Qk_hat + SDR_term_k_cum)
            # mean(data$dat.sVar[t==0, Qk_hat])
            # mean(data$dat.sVar[t==0, Qk_hat] + SDR_term_k)
            # mean(data$dat.sVar[t==0, Qk_hat] + SDR_term_k_cum)
            # # mean(Qk_hat + SDR_term_k_cum)
            # mean(SDR_term_k_cum)
            # mean(SDR_term_k)
            # i <- 2
            # private$PsAsW.models[[i]]$transform_Q_k(data = data,
            #                                           k_i = k_i,
            #                                           max_Qk_idx = max_Qk_idx,
            #                                           ...)
            SDR_term_k <- private$PsAsW.models[[i]]$transform_Q_k(data = data,
                                                    k_i = k_i,
                                                    i = i,
                                                    max_Qk_idx = max_Qk_idx,
                                                    SDR_term_k_cum = SDR_term_k_cum,
                                                    # QModel_h_k = QModel_h_k,
                                                    # QModel_Qkplus1 = QModel_Qkplus1,
                                                    ...)
            SDR_term_k_cum <- SDR_term_k_cum + SDR_term_k
            # colMeans(cbind(SDR_term_k_cum, SDR_term_k, SDR_term_k_cum + SDR_term_k))
            # Qk_hat <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
            # print("mean(SDR_term_k_cum): "); print(mean(SDR_term_k_cum))
            # print("\nMean estimand of E(Y_d), t = 1: ")
            # print(mean(data$dat.sVar[t==1, Qk_hat]))
            # print("\nMean estimand of E(Y_d), t = 0: ")
            # print(mean(data$dat.sVar[t==0, Qk_hat]))
            # print("\nMean Qkplus1 of E(Y_d), t = 0: ")
            # print(mean(data$dat.sVar[t==0, Qkplus1]))
            # print("DR transform estimand of E(Y_d), t = 1: ")
            # print(mean(data$dat.sVar[t==1, Qk_hat] +  SDR_term_k_cum))
            # print("DR transform estimand of E(Y_d), t = 0: ")
            # print(mean(data$dat.sVar[t==0, Qk_hat] +  SDR_term_k_cum))

            # mean(data$dat.sVar[t==0, Qk_hat])
            # glm.fit <- glm(Qk_hat ~ 1, data = data$dat.sVar[t==0, "Qk_hat"])

            # browser()
          # }
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
    transform_Q_k = function(data, k_i, i, max_Qk_idx, SDR_term_k_cum = 0, ...) {
      ## use only the observations that participated in fitting of the initial Q_{k'} (current time-point is k')
      use_subset_idx <- self$idx_used_to_fit_initQ
      # if (gvars$verbose == 2) {
        cat("...running SDR targeting loop...\n")
        cat("Total number of time-points to consider: " %+% max_Qk_idx, "\n")
        cat("Current (S)DR k index = " %+% self$all_Qregs_indx[k_i], "\n")
        cat("Current (S)DR k' index = " %+% self$all_Qregs_indx[i], "\n")
      # }

      data$IPwts_by_regimen[use_subset_idx,][1:100, ]
      data$IPwts_by_regimen[self$subset_idx, ][1:100, ]

      data$dat.sVar[use_subset_idx,][1:100, ]
      data$dat.sVar[self$subset_idx,][1:100, ]

      ## 1. Weights: defined new column of cumulative weights where cumulative product starts at t = Qk_idx (k), rather than t = 0:
      # wts <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      wts <- data$IPwts_by_regimen[self$subset_idx, "cum.IPAW", with = FALSE][[1]]
      ## 2a. Outcome: **TARGETED** prediction of the previous step k'+1.
      # Qkplus1.protected <- data$dat.sVar[use_subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      # Qkplus1.protected <- data$dat.sVar[self$subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      Qkplus1.protected <- data$dat.sVar[self$subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      ## 3. Initial model fits:
      # Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]
      Qk_hat <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
      ## 3. Initial model fit that includes extrapolated predictions (newsly censored and non-rule followers):
      # Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]

      ## The EIC estimates and the cumulative transform sum:
      ## REPLACE Qk_hat with Qk_hat_transform (targeted version)
      SDR_term_k <- wts * (Qkplus1.protected - Qk_hat)
      print("mean(SDR_term_k)"); print(mean(SDR_term_k))

      wts_fitted <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      Qkplus1.protected_fitted <- data$dat.sVar[use_subset_idx, "Qkplus1.protected", with = FALSE][[1]]
      Qk_hat_fitted <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]
      SDR_term_k_fitted <- wts_fitted * (Qkplus1.protected_fitted - Qk_hat_fitted)
      print("mean(SDR_term_k_fitted)"); print(mean(SDR_term_k_fitted))
      # Qk_hat_fitted_2 <- Qk_hat_fitted + SDR_term_k_fitted
      # mean(wts_fitted * (Qkplus1.protected_fitted - Qk_hat_fitted_2))
      # mean(Qk_hat)
      # mean(Qkplus1.protected)
      # mean(SDR_term_k)
      # browser()

      # data$dat.sVar[use_subset_idx, ("EIC_i_t") := SDR_term_k]
      # data$dat.sVar[use_subset_idx, EIC_i_t_sum]
      # data$dat.sVar[use_subset_idx,]
      # browser()
      # data$dat.sVar[use_subset_idx, ("EIC_i_t_sum") := EIC_i_t_sum + SDR_term_k_fitted]
      data$dat.sVar[self$subset_idx, ("EIC_i_t_sum") := EIC_i_t_sum + SDR_term_k]
      # data$dat.sVar[self$subset_idx, EIC_i_t_sum := 0]
      # mean(data$dat.sVar[self$subset_idx-1,Qkplus1.protected])
      # mean(data$dat.sVar[self$subset_idx-1, Qk_hat])

      ## Signifies we are in a final (S)DR transform loop (final / last regression for t=0 has been reached).
      ## No more fits of the initials are going to be performed, just evaluate the EIC
      if (k_i == max_Qk_idx) {
        data$dat.sVar[use_subset_idx, ("EIC_i_t") := SDR_term_k_fitted]
        # data$dat.sVar[self$subset_idx, ("EIC_i_t") := SDR_term_k]
      }
      ## Signifies we reached a final iteration of the final (S)DR transform loop. No more initial fits are going to be performed.
      if (i == max_Qk_idx) {

      }

      ## Reached the last iter of current SDR loop.
      ## This is where we actually define (S)DR transform Q^*[k].
      ## This is our current model update and it will be used as the outcome for the next regression.
      if (i == k_i) {
        # Qk_hat <- Qk_hat + data$dat.sVar[use_subset_idx, EIC_i_t_sum]
        # Qk_hat <- Qk_hat + SDR_term_k
        Gamma_DR <- Qk_hat + data$dat.sVar[self$subset_idx, EIC_i_t_sum]
        # Gamma_DR <- Qk_hat_fitted + data$dat.sVar[use_subset_idx, EIC_i_t_sum]

        # data$dat.sVar[use_subset_idx, "Qk_hat" := Gamma_DR]
        # data$dat.sVar[self$subset_idx, "Qk_hat" := Gamma_DR]

        # mean(Gamma_DR)
        # data$dat.sVar[self$subset_idx, ]

      }

      # data$dat.sVar[use_subset_idx,]
      # data$dat.sVar[self$subset_idx,]
      # data$dat.sVar[use_subset_idx - 1,]
      # data$dat.sVar[self$subset_idx - 1,]
      # browser()

      # Over-write the old predictions with new model updates as Qk_hat[k'] in row [k']:
      # data$dat.sVar[self$subset_idx, "Qk_hat" := Gamma_DR]
      # Gamma_DR-data$dat.sVar[self$subset_idx, Qk_hat]
      # print("MSE of previous Qk_hat vs. upated Qk_hat: " %+% mean((Gamma_DR-Qk_hat_all)^2))

      if ((self$Qreg_counter > 1) && (i == k_i)) {
        # if (gvars$verbose)
        # Qk_hat_transform <- Qk_hat + data$dat.sVar[use_subset_idx, EIC_i_t_sum]
        cat("updating the outcomes for next regression with DR tranformed outcomes (only among those who were at risk and following treatment)")
        # data$dat.sVar[use_subset_idx - 1, "Qkplus1.protected" := Qk_hat]
        # data$dat.sVar[self$subset_idx - 1, "Qkplus1.protected" := data$dat.sVar[self$subset_idx, "Qk_hat"]]
        data$dat.sVar[self$subset_idx - 1, "Qkplus1.protected" := Qkplus1]
        data$dat.sVar[self$subset_idx - 1, "Qkplus1" := Gamma_DR]
        # data$dat.sVar[self$subset_idx - 1, "Qkplus1" := data$dat.sVar[self$subset_idx, Qk_hat]]
        # data$dat.sVar[use_subset_idx - 1, "Qkplus1" := Gamma_DR]

      } else if (self$Qreg_counter > 1) {
        ## Q: Is this right????
        ## I.e., is it possible to have subjects that previously contributed to rolling sum, but have no contributed
        ## at this time-point?
        # data$dat.sVar[self$subset_idx - 1, ("Qk_hat") := Qk_hat_transform]
        data$dat.sVar[self$subset_idx - 1, ("EIC_i_t_sum") := data$dat.sVar[self$subset_idx, EIC_i_t_sum]]

      } else {
        # browser()
        data$dat.sVar[self$subset_idx, ("Qk_hat") := Gamma_DR]
        # private$probAeqa <- data$dat.sVar[self$subset_idx, Qk_hat]
        private$probAeqa <- Gamma_DR
        # private$probAeqa <- Gamma_DR
        # data$dat.sVar[]
        # mean(Gamma_DR)
      }

      if (i == k_i) {
        print("mean(Gamma_DR)"); print(mean(Gamma_DR))
        return(Gamma_DR)
      } else {
        return(SDR_term_k)
      }

      # return(SDR_term_k)
      # invisible(self)
    }
  )
)