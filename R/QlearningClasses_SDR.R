
##**********************************************************************
################ TO DO FOR SDR ################
##**********************************************************************
##
## * 1. Remove manual SDR update fitting / prediction and replace with general fit / predict for any learner.
##      Make it consistent with fit / predict in the rest of the package.
##      Allow passing arbitrary learners for the targeted update via model syntax.
##
## * 2. Stratification by A[k-1]? Need to sort out how to do this correctly for all targeting steps.
##
## * 3. Do we need to use offsets when extrapolating predictions of Q.star.k for new observations
##      (newly censored and new non-follows at k')?
##
## * (DONE) Safe targeting of E[Y_d] without over-writing E[Y_d|W].
##      When k=-1 (after the last run of the epsilon targeting regressions).
##      Last loop can be currently identified as (Qk_idx == max_Qk_idx).
##      This final targeting step will do the usual TMLE (intercept-only logistic regression updates).
##      NEED TO FIGURE OUT HOW ACCOMPLISHE THIS AND WHERE TO STORE THE PREDICTIONS
##      (without over-overwriting the targeted E[Y_d|W] from previous loop)
##
## * (DONE) This cannot conflict with predictions for E[Y_d]. I.e.,
##      Should be able to make prediction from E(Y_d|W') and obtain E(Y_d) from same run of the algorithm.
##      Look at regular GCOMP for analogy on how best to set this up.
##
## * (DONE) Set up for doing predictions for targeted functional E(Y_d|W') for any W'. This needs to be entirely automated.
## * (DONE) Finish the weights evaluation, see below (need to pass appropriate args to the class to make the function call to weights).
## * (DONE) Still need to figure out how to pass covariates from the current G-COMP fit of Q_k, the covariates need to be from k-1.
## * (DONE) Add newly censored obs and newly non-followers to targeted predictions.
## **********************************************************************


## ---------------------------------------------------------------------
## R6 class for fitting SDR procedure
## Inherits from \code{GenericModel}.
## ---------------------------------------------------------------------
SDRModel <- R6Class(classname = "SDRModel",
  inherit = GenericModel,
  portable = TRUE,
  class = TRUE,
  public = list(

    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      # serial loop over all regressions in PsAsW.models:
      max_Qk_idx <- length(private$PsAsW.models)
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$fit(data = data, ...)

        ## The last targeting step (k_i==max_Qk_idx) should be the usual TMLE and it needs to be done differently
        ## (univariate logistic model updates)
        for (i in (1:k_i)) {
          ## All the targeting is for one functional with reg index (Qk_idx + 1) (since time-points and loops are reversed)
          ## Thus, all the targeting steps in this loop are functions of the same covariate space, located in row shift:
          ## Hk_row_offset = kprime_idx - (Qk_idx + 1)

          ## Qk_idx = 1, kprime_idx = 1 -> Hk_row_offset = -1

          ## Qk_idx = 2, kprime_idx = 1 -> Hk_row_offset = -2
          ## Qk_idx = 2, kprime_idx = 2 -> Hk_row_offset = -1

          ## Qk_idx = 3, kprime_idx = 1 -> Hk_row_offset = -3
          ## Qk_idx = 3, kprime_idx = 2 -> Hk_row_offset = -2
          ## Qk_idx = 3, kprime_idx = 3 -> Hk_row_offset = -1

          ## Qk_idx = 3, kprime_idx = 1 -> Hk_row_offset = -3
          ## Qk_idx = 3, kprime_idx = 2 -> Hk_row_offset = -2
          ## Qk_idx = 3, kprime_idx = 3 -> Hk_row_offset = -1

          private$PsAsW.models[[i]]$eval_weights_k(data = data, ...)

          private$PsAsW.models[[i]]$fit_epsilon_Q_k(data = data,
                                                    kprime_idx = i,
                                                    Qk_idx = k_i,
                                                    max_Qk_idx = max_Qk_idx, ...)
        }
      }
      invisible(self)
    }
  )
)

## ---------------------------------------------------------------------
## R6 Class for Sequentially Double Robustness Targeting Procedure
## Internal implementation of Q-learning functionality.
## Inherits from \code{QlearnModel} R6 Class.
## ---------------------------------------------------------------------
SDRQlearnModel  <- R6Class(classname = "SDRQlearnModel",
  inherit = QlearnModel,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    regimen_names = character(), # for future pooling across regimens
    wts_k = NULL,
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,
    Qreg_counter = integer(), # Counter for the current sequential Q-regression (min is at 1)
    t_period = integer(),

    ## **********************************************************************
    ## Weights for SDR
    ## **********************************************************************
    ## *THE VALUE OF the current time-point is saved in self$t_period.
    ## * When this is passed as tmin to getIPWeights, the weights are evaluated correctly for current k.
    ## * That is, the product of the cumulative weights will be taken starting with k=t all the way up to end of follow-up.
    ## * The way the weights need to be evaluated for each k is exactly like in getIPWeights(),
    ##    ****** except that all the weights for time-points < tmin are set to constant 1.
    eval_weights_k = function(data, ...) {
      call_list <- data$IPWeights_info
      call_list[["tmin"]] <- self$t_period
      call_list[["OData"]] <- data
      if (gvars$verbose == 2) message("...evaluating IPWeights for SDR...")
      IPWeights <- do.call("getIPWeights", call_list)
      data$IPwts_by_regimen <- IPWeights
      return(invisible(IPWeights))
    },

    ## This update is for current time-point k' aimed at targeting initial Q[k] fit.
    ## Outcome: The same outcomes that were used to fit this initial Q[k']
    ## Predictors: All predictors used for fitting Q[k-1]
    ## Weights: the product of g's from t=k to current k'
    ## Offset: qlogis(Qk_hat) - current (initial) Q prediction (at k')
    fit_epsilon_Q_k = function(data, kprime_idx, Qk_idx, max_Qk_idx, ...) {
      if (self$all_Qregs_indx[kprime_idx] != self$Qreg_counter) stop("something terrible has happened")

      ## All targeting is for one functional with loop index (Qk_idx + 1) (since loops are reverse of time-points)
      ## Thus, all the targeting steps in this loop are functions of the same covariate space, located in row shift:
      Hk_row_offset <- kprime_idx - (Qk_idx + 1)

      ## use only the observations that participated in fitting of the initial Q_{k'} (current time-point is k')
      use_subset_idx <- self$idx_used_to_fit_initQ

      if (gvars$verbose == 2) {
        cat("...running SDR targeting loop...\n")
        cat("Total number of time-points: " %+% max_Qk_idx, "\n")
        cat("Updating Q_k at index: " %+% self$all_Qregs_indx[Qk_idx], "\n")
        cat("Current targeting step for Q_k' is at index: " %+% self$Qreg_counter, "; time-point: " %+% self$t_period, "\n")

        ## above is the same as self$all_Qregs_indx[kprime_idx]
        cat("Current k' (kprime) idx = " %+% kprime_idx, "\n")
        cat("Currently targeting covariate space in time-point (kprime) idx: " %+% kprime_idx, "\n")

        cat("The shift row for targeted covariate space: " %+% Hk_row_offset, "\n")
        cat("Targeting covariate space for time-point: " %+% (self$t_period+Hk_row_offset), "\n")
        cat("length(use_subset_idx): ", length(use_subset_idx), "length(self$subset_idx): ", length(self$subset_idx), "\n")
      }

      ## 1. Weights: defined new column of cumulative weights where cumulative product starts at t = Qk_idx (k), rather than t = 0:
      wts <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]

      ## 2. Outcome: **TARGETED** prediction of the previous step k'+1.
      ##    These were targeted towards estimation of the same Q_k by the previous call to this function.
      ##    If this is the first call to this function, then these must be the initial outcomes. That is, Qkplus1 is INITIALIZED TO BE ALL Y AT FIRST.
      ##    This is used as the outcome for the current regression update.
      Qkplus1 <- data$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]

      ## 3. Offset:
      ##    Our offset is the Q[k'] fit that was previously targeted towards estimation of Q[k'] (k' is current time point).
      ##    This estimand has also been saved in row k'-1 by previous initial G-COMP for k'.
      ##    We now re-target it towards estimation of Q[k] (for any k < k').
      ##    In standard LTMLE this role is played by the initial GCOMP prediction for Q[k'].
      ##    If this is the first call to this function for k' then this is exactly the initial GCOMP prediction.
      ##    Note that after we generate an update Qk_hat_star, we ****over-write**** current Qk_hat with its updated version
      Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]

      ## 4A. The model update. Univariate logistic regression (TMLE)
      if (Qk_idx == max_Qk_idx) {
        # browser()
        if (gvars$verbose) cat("Last targeting step for E(Y_d) with intercept only TMLE updates\n")
        # Qk_hat_star_all <- intercept.update(data,
        #                                     Qkplus1 = Qkplus1,
        #                                     Qk_hat = Qk_hat,
        #                                     wts = wts,
        #                                     lower_bound_zero_Q = self$lower_bound_zero_Q,
        #                                     skip_update_zero_Q = self$skip_update_zero_Q)
        TMLE.fit <- tmle.update(Qkplus1 = Qkplus1,
                                Qk_hat = Qk_hat,
                                IPWts = wts,
                                lower_bound_zero_Q = self$lower_bound_zero_Q,
                                skip_update_zero_Q = self$skip_update_zero_Q)
        TMLE_intercept <- TMLE.fit$TMLE_intercept
        if (!is.na(TMLE_intercept) && !is.nan(TMLE_intercept)) {
          update.Qstar.coef <- TMLE_intercept
        } else {
          update.Qstar.coef <- 0
        }
        ## Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
        Qk_hat_star_all <- plogis(qlogis(Qk_hat_all) + update.Qstar.coef)

      ## 4B. The model update. Infinite dimensional epsilon (SDR)
      } else {
        # Qk_hat_star_all <- SDR.update()
        ## 1. Grabbing covariates for H(k-1), where k-1 is determined by the row index offset Hk_row_offset
        ##    somehow we need to select from the same set of covariates only those people who are currently at risk
        ##    => WE NEED A WAY OF GRABING THE stratification (init_idx) THAT WERE USED (OR WILL BE USED) FOR PREDICTING THE INITAL FOR TIME-POINT k-1.
        ##    Then we can intersect those with current idx
        obs_dat <- data$dat.sVar[use_subset_idx + Hk_row_offset, self$predvars, with = FALSE]
        # data$dat.sVar[t == 10, ]; data$dat.sVar[t == 9, ]

        ## 2. Fitting regression: Qkplus1 ~ offset(qlogis(Qk_hat)) + H[k-1] and weights 'wts'
        # require('xgboost')
        # params <- list("objective" = "reg:logistic", "booster" = "gbtree", "nthread" = 1, "max_delta_step" = 10)
        params <- self$reg$SDR_model
        cat("running SDR w/ following params: \n "); str(params)
        # obs_dat[, CVD := as.numeric(CVD)]
        xgb_dat <- xgboost::xgb.DMatrix(as.matrix(obs_dat), label = Qkplus1)
        xgboost::setinfo(xgb_dat, "base_margin", qlogis(Qk_hat))
        xgboost::setinfo(xgb_dat, "weight", wts)

        nrounds <- params[["nrounds"]]
        params[["nrounds"]] <- NULL
        # browser()
        if (is.null(nrounds)) {
          cat("...running cv to figure out best nrounds for epsilon target...\n")
          mfitcv <- xgboost::xgb.cv(params = params, data = xgb_dat, nrounds = 100, nfold = 5, early_stopping_rounds = 10)
          nrounds <- mfitcv$best_iteration
          cat("...best nrounds: ", nrounds, "\n")
        }
        mfit <- xgboost::xgb.train(params = params, data = xgb_dat, nrounds = nrounds)

        # require('gam')
        # mfit <- gam.fit(
        #   x=as.matrix(data$dat.sVar[use_subset_idx,self$predvars]),
        #   y = Qstarkprime,
        #   family=binomial(),
        #   offset = offset(qlogis(Qk_hat)),
        #   weights=wts)
        # pred1.star = as.numeric(predict(mfit,type='response'))

        ## 3. Predicting for all newly censored and new non-follows at k'
        ##    **** Q: DO WE EVEN NEED TO USE NEW OFFSETS WHEN WE ARE EXTRAPOLATING THE MODEL UPDATE????
        ##    **** Q: HOW TO GENERALIZE THIS TO MAKE PREDICTIONS of E[Y_d|W'] for any new W'?
        ##    **** NEED TO BE ABLE TO MAKE UPDATED MODEL PREDICTIONS FOR ANY NEW DATASET BASED ON THIS MODEL UPDATE
        ##    THE wts used no longer play any role, but the offsets (Qk_hat) needs to be re-evaluated for new observations?
        ##    The predictors used are still based on the covariate space at time point k-1.
        pred_dat <- data$dat.sVar[self$subset_idx + Hk_row_offset, self$predvars, with = FALSE]
        # pred_dat[, CVD := as.numeric(CVD)]
        xgb_dat <- xgboost::xgb.DMatrix(as.matrix(pred_dat))
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
        xgboost::setinfo(xgb_dat, "base_margin", qlogis(Qk_hat_all))

        ## 4. Update the model predictions (Qk_hat) for initial Q[k'] from GCOMP at time-point k'.
        ##    Based on TMLE update, the predictions now include ALL obs that are newly censored
        ##    and just stopped following the rule at current k'.
        Qk_hat_star_all <- predict(mfit, xgb_dat)
      }

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