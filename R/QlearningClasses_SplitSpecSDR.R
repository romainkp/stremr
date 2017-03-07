## ---------------------------------------------------------------------
## R6 Class for Split-Specific Sequentially Double Robustness Targeting Procedure
## Internal implementation of Q-learning functionality.
## Inherits from \code{QlearnModel} R6 Class.
## ---------------------------------------------------------------------
SplitCVSDRQlearnModel  <- R6Class(classname = "SplitCVSDRQlearnModel",
  inherit = SDRQlearnModel,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(

    ## Split-Specific SL updater
    ## This update is for current time-point k' aimed at targeting initial Q[k] fit.
    ## Outcome: The same outcomes that were used to fit this initial Q[k']
    ## Predictors: All predictors used for fitting Q[k-1]
    ## Weights: the product of g's from t=k to current k'
    ## Offset: qlogis(Qk_hat) - current (initial) Q prediction (at k')
    fit_epsilon_Q_k = function(data, kprime_idx, Qk_idx, max_Qk_idx, QModel_h_k, ...) {
      if (self$all_Qregs_indx[kprime_idx] != self$Qreg_counter) stop("something terrible has happened")
      ## All targeting is for one functional with loop index (Qk_idx + 1) (since loops are reverse of time-points)
      ## Thus, all the targeting steps in this loop are functions of the same covariate space, located in row shift:
      Hk_row_offset <- kprime_idx - (Qk_idx + 1)
      ## use only the observations that participated in fitting of the initial Q_{k'} (current time-point is k')
      use_subset_idx <- self$idx_used_to_fit_initQ

      if (gvars$verbose == 2) {
        cat("...running SDR targeting loop...\n")
        cat("Total number of time-points to consider: " %+% max_Qk_idx, "\n")
        cat("Current SDR k' (kprime) index = " %+% self$all_Qregs_indx[kprime_idx], "\n")
        cat("Targeting the same Q.kplus1 for t index = " %+% self$all_Qregs_indx[Qk_idx], "\n")
        cat("Using the covariate space at t index = " %+% (self$all_Qregs_indx[kprime_idx]+Hk_row_offset), "\n")
      }

      ## 1. Weights: defined new column of cumulative weights where cumulative product starts at t = Qk_idx (k), rather than t = 0:
      wts <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      ## 2. Outcome: **TARGETED** prediction of the previous step k'+1.
      Qkplus1 <- data$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]
      ## 3. Offset:
      Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]

      ## 4A. The model update. Univariate logistic regression (TMLE)
      if (Qk_idx == max_Qk_idx) {
        if (gvars$verbose) cat("Last targeting step for E(Y_d) with intercept only TMLE updates\n")
        ## Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
        X <- as.matrix(qlogis(Qk_hat))
        newX <- as.matrix(qlogis(Qk_hat_all))
        colnames(X) <- colnames(newX) <- "offset"
        # TMLE.fit <- SDR.updater.speedglmTMLE(Y = Qkplus1, X = X, newX = newX, obsWeights = wts)
        TMLE.fit <- SDR.updater.glmTMLE(Y = Qkplus1, X = X, newX = newX, obsWeights = wts)
        Qk_hat_star_all <- TMLE.fit[["pred"]]
        EIC_i_t_calc <- wts * (Qkplus1 - Qk_hat)
        data$dat.sVar[use_subset_idx, ("EIC_i_t") := EIC_i_t_calc]

      ## 4B. The model update. Infinite dimensional epsilon (SDR)
      } else {
        ## Fitting SDR update with split-spec SL: Qkplus1 ~ offset(qlogis(Qk_hat)) + H[k-1] and weights 'wts'
        ## Grabbing covariates for H(k-1), where k-1 is determined by the row index offset Hk_row_offset
        ## Predict for all newly censored and new non-follows at k'
        ##    GENERALIZED TO MAKE PREDICTIONS of E[Y_d|W'] for any new W'  BASED ON THIS MODEL UPDATE
        ##    THE wts used no longer play any role, but the offsets (Qk_hat) needs to be re-evaluated for new observations
        ##    The predictors used are still based on the covariate space at time point k-1.
        ##    Update the model predictions (Qk_hat) for initial Q[k'] from GCOMP at time-point k'.
        ##    Based on TMLE update, the predictions now include ALL obs that are newly censored
        ##    and just stopped following the rule at current k'.

        ##  Grab covariates for H(k-1), where k-1 is determined by the row index offset Hk_row_offset
        epsilon_predvars <- QModel_h_k$predvars
        obs_dat <- data$dat.sVar[use_subset_idx + Hk_row_offset, epsilon_predvars, with = FALSE]
        obs_dat[, ("offset") := qlogis(Qk_hat)]

        pred_dat <- data$dat.sVar[self$subset_idx + Hk_row_offset, epsilon_predvars, with = FALSE]
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
        pred_dat[, ("offset") := qlogis(Qk_hat_all)]


        ## ------------------------------------------------------------------------------
        ## ****** MOVE THIS OUTSIDE (NEED TO PASS folds OBJECT ONLY) ******
        ## ------------------------------------------------------------------------------
        make_kfold_from_column <- function(data, id = ".id", fold_column = "fold") {
          n <- nrow(data)
          folds <- data[[fold_column]]
          k <- length(unique(folds))

          idx <- seq_len(n)
          fold_idx <- split(idx, folds)

          fold <- function(v, test) {
              origami:::make_fold(v, setdiff(idx, test), test)
          }
          purrr::map2((1:k), fold_idx, fold)
        }
        # make_kfold_from_column
        folds <- make_kfold_from_column(data$dat.sVar[use_subset_idx + Hk_row_offset], id = ".id", fold_column = "fold_ID")

        # cv_split_preds <- function(fold, data, fits_Qk, use_full = FALSE) {
        #     browser()
        #     ## These will be automatically defined in the calling frame of this function
        #     ## when the cross-validator that calls cv_split_preds()
        #     v <- origami::fold_index()
        #     train_idx <- origami::training()
        #     valid_idx <- origami::validation()

        #     if ((all.equal(train_idx, valid_idx) == TRUE) || use_full) {
        #       # we're in final resubstitution call, so let's use full Q and g fits
        #       splitQk_fit <- fits_Qk$fullFit
        #     } else {
        #       # split-specific Super Learners
        #       splitQk_fit <- fits_Qk$foldFits[[v]]
        #     }
        #     ## split-specific predictions for new data
        #     ## this may include new observations, e.g., extrapolating for newly censored
        #     new_data <- as.matrix(data)
        #     QAW <- predict(splitQ_fit, newdata = new_data)[["pred"]]
        #     # new_data[, nodes$Anode] <- 0
        #     # Q0W <- predict(splitQ_fit, newdata = new_data)$pred
        #     # new_data[, nodes$Anode] <- 1
        #     # Q1W <- predict(splitQ_fit, newdata = new_data)$pred
        #     # pA1 <- predict(splitg_fit, new_data)$pred

        #     # # split specific blip, class, and weights
        #     # A <- data[, nodes$Anode]
        #     # Y <- data[, nodes$Ynode]
        #     # D1 <- (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W
        #     # if (maximize) {
        #     #     Z <- as.numeric(D1 > 0)
        #     # } else {
        #     #     Z <- as.numeric(D1 < 0)
        #     # }
        #     # K <- as.vector(abs(D1))  #D1 is a matrix somehow

        #     # browser()

        #     list(QAW = QAW, pA1 = pA1)
        # }


        ## PASS THE SPLIT-SPEC PREDS FROM THE PREVIOUS RUN
        ## NEED TO EXTRACT SPLIT-SPEC Y and offset
        SL.library <- c("SDR.updater.NULL", "SDR.updater.glmTMLE", "SDR.updater.glm", "SDR.updater.xgb")
        # , "SDR.updater.speedglmTMLE"
        library("abind")
        # browser()
        SDR_SL_fit <- origami::origami_SuperLearner(folds = folds,
                                                    Y = Qkplus1,
                                                    X = as.matrix(obs_dat),
                                                    family = quasibinomial(),
                                                    obsWeights = wts,
                                                    SL.library = SL.library,
                                                    params = self$reg$SDR_model)
        print("SDR_SL_fit: "); print(SDR_SL_fit)
        ## SuperLearner final prediction based on models fit on all data
        Qk_hat_star_all <- as.numeric(predict(SDR_SL_fit, as.matrix(pred_dat))[["pred"]])
        # mfit <- SDR.updater.xgb(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts, params = self$reg$SDR_model)
        # # mfit <- SDR.updater.glm(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # # mfit <- SDR.updater.TMLE(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # # mfit <- SDR.updater.NULL(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # Qk_hat_star_all <- mfit[["pred"]]


        # ## Split-Specific predictions from the SuperLearner trained on fold v, for new data (extrapolating to new obs)
        # v <- 1
        # splitQk_fit <- SDR_SL_fit$foldFits[[v]]
        # split_preds_v <- as.numeric(predict(splitQk_fit, newdata = as.matrix(pred_dat))[["pred"]])
        # split_preds <- origami::cross_validate(cv_split_preds, folds, pred_dat, fits, .combine = FALSE)
        # browser()
        # # names(SDR_SL_fit)
        # # SDR_SL_fit[["foldFits"]][[1]]
        # # SLpred <- predict(SDR_SL_fit, as.matrix(pred_dat))



      }

      print("MSE of previous Qk_hat vs. upated Qk_hat: " %+% mean((Qk_hat_star_all-Qk_hat_all)^2))

      # Over-write the old predictions with new model updates as Qk_hat[k'] in row [k']:
      data$dat.sVar[self$subset_idx, "Qk_hat" := Qk_hat_star_all]

      ## ************ NOTE ************
      ##  The shift below will be wrong when we get to the very last iteration of the very last update (LTMLE)
      ##  There will be nowhere to write a new update, since we will be at a row t==0.
      ##  This happens when: (Qk_idx == max_Qk_idx) && (Qk_idx == kprime_idx)
      ## ************
      ## Set the outcome for the next Q-regression: put the updated Q[k'] in (k'-1)
      ## (kprime_idx == Qk_idx) means that we have reached the final update for current Q[k] that we were targeting,
      ## Save this update in row [k'-1] = [k-1], that's where it will be picked up by next loop (or initial est step for Q[k-1]).
      if ((Qk_idx == max_Qk_idx) && (Qk_idx == kprime_idx)) {
        if (gvars$verbose) cat("reached the last targeting iteration of the very last initial regression. Saving the final targeted prediction for E[Y_d]")
        private$probAeqa <- Qk_hat_star_all
      } else {
        data$dat.sVar[(self$subset_idx - 1), "Qkplus1" := Qk_hat_star_all]
      }
      invisible(self)
    }
  )
)