
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

## **********************************************************************
## Weights for SDR
## **********************************************************************
## *THE VALUE OF the current time-point is saved in self$t_period.
## * When this is passed as ignore_tmin to getIPWeights, the weights are evaluated correctly for current k.
## * That is, the product of the cumulative weights will be taken starting with k=t all the way up to end of follow-up.
## * The way the weights need to be evaluated for each k is exactly like in getIPWeights(),
##    ****** except that all the weights for time-points < ignore_tmin are set to constant 1.
eval_weights_k = function(data, ignore_tmin = NULL, ignore_max = NULL, reverse_wt_prod = FALSE, ...) {
  call_list <- data$IPWeights_info

  if (!is.null(ignore_tmin)) {
    call_list[["ignore_tmin"]] <- ignore_tmin
  }
  # else {
  #   call_list[["ignore_tmin"]] <- self$t_period
  # }

  if (!is.null(ignore_max)) {
    call_list[["ignore_max"]] <- ignore_max
  }

  call_list[["reverse_wt_prod"]] <- reverse_wt_prod

  call_list[["OData"]] <- data
  if (gvars$verbose == 2) message("...evaluating IPWeights for SDR...")

  IPWeights <- do.call("getIPWeights", call_list)
  data$IPwts_by_regimen <- IPWeights
  return(invisible(IPWeights))
}

## ---------------------------------------------------------------------
## R6 class for fitting SDR procedure
## Inherits from \code{ModelGeneric}.
## ---------------------------------------------------------------------
SDRModel <- R6Class(classname = "SDRModel",
  inherit = ModelGeneric,
  portable = TRUE,
  class = TRUE,
  public = list(

    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      # serial loop over all regressions in PsAsW.models:
      max_Qk_idx <- length(private$PsAsW.models)
      for (k_i in seq_along(private$PsAsW.models)) {

        private$PsAsW.models[[k_i]]$fit(data = data, ...)

        ## The inner loop targets Q.kplus1, but based on the same time-point (covariate space) for object (k_i+1).
        ## If (k_i+1) doesn't exist, it means that we have reached the initial regression E(Q.kplus1|A,W).
        ## This requires the last targeting step, which is different from the rest.
        ## The last targeting step (when k_i==max_Qk_idx) should be the usual TMLE update (univariate logistic model updates).


        ## 1. Previous SDR updating scheme, loop the updates over different Q[i], wrt to the same covariate space of h[k_i-1] (object private$PsAsW.models[[k_i + 1]])
        # if (k_i == max_Qk_idx) {
        #   QModel_h_k <- NULL
        # } else {
        #   QModel_h_k <- private$PsAsW.models[[k_i + 1]]
        # }
        # for (i in (1:k_i)) {

        ## 2. New SDR targeting scheme, loop by updating the same Q[k+i], wrt to different covariate spaces of h[i]
        for (i in (k_i:max_Qk_idx)) {
          if (i == max_Qk_idx) {
            QModel_h_k <- NULL
          } else {
            QModel_h_k <- private$PsAsW.models[[i + 1]]
          }

          if (k_i==1L) {
            QModel_Qkplus1 <- NULL
          } else {
            QModel_Qkplus1 <- private$PsAsW.models[[k_i-1]]
          }
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

          ## 1. Previous SDR targeting scheme:
          ## evaluate the weights for this targeting step:
          # #### private$PsAsW.models[[i]]$eval_weights_k(data = data, , ignore_tmin = private$PsAsW.models[[i]]$t_period, ...)
          # wts <- eval_weights_k(data = data, ignore_tmin = private$PsAsW.models[[i]]$t_period, ...)
          # private$PsAsW.models[[i]]$fit_epsilon_Q_k(data = data,
          #                                           kprime_idx = i,
          #                                           Qk_idx = k_i,
          #                                           max_Qk_idx = max_Qk_idx,
          #                                           QModel_h_k = QModel_h_k,
          #                                           QModel_Qkplus1 = QModel_Qkplus1,
          #                                           ...)


          ## 2. New SDR targeting scheme. Always updating the very same Q[k_i] we just fit as initial
          ## ******** Q: NEED TO VERIFY ignore_tmin is ACTUALLY SET-UP CORRECTLY ********
          # private$PsAsW.models[[k_i]]$eval_weights_k(data = data, ignore_tmin = private$PsAsW.models[[i]]$t_period, ...)
          wts <- eval_weights_k(data = data, ignore_tmin = private$PsAsW.models[[i]]$t_period, ...)
          ## ********

          private$PsAsW.models[[k_i]]$fit_epsilon_Q_k(data = data,
                                                      kprime_idx = k_i,
                                                      Qk_idx = i,
                                                      max_Qk_idx = max_Qk_idx,
                                                      QModel_h_k = QModel_h_k,
                                                      QModel_Qkplus1 = QModel_Qkplus1,
                                                      ...)
        }
      }
      invisible(self)
    }
  )
)

## ---------------------------------------------------------------------
## R6 Class for Sequentially Double Robustness Targeting Procedure
## Internal implementation of Q-learning functionality.
## Inherits from \code{ModelQlearn} R6 Class.
## ---------------------------------------------------------------------
SDRModelQlearn  <- R6Class(classname = "SDRModelQlearn",
  inherit = ModelQlearn,
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

    ## This update is for current time-point k' aimed at targeting initial Q[k] fit.
    ## Outcome: The same outcomes that were used to fit this initial Q[k']
    ## Predictors: All predictors used for fitting Q[k-1]
    ## Weights: the product of g's from t=k to current k'
    ## Offset: qlogis(Qk_hat) - current (initial) Q prediction (at k')
    fit_epsilon_Q_k = function(data, kprime_idx, Qk_idx, max_Qk_idx, QModel_h_k, QModel_Qkplus1, ...) {
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
      ##    These were targeted towards estimation of the same Q_k by the previous call to this function.
      ##    If this is the first call to this function, then these must be the initial predictions. Note that Qkplus1 is INITIALIZED TO BE ALL Y AT FIRST.
      ##    This is used as the outcome for the current regression update.
      Qkplus1 <- data$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]

      ## 3. Offset:
      ##    Our offset is the Q[k'] fit that was previously targeted towards estimation of Q[k'] (k' is the current time point).
      ##    This estimand has also been saved in row k'-1 by previous initial G-COMP for k'.
      ##    We now re-target it towards estimation of Q[k] (for any k < k').
      ##    In standard LTMLE this role is played by the initial GCOMP prediction for Q[k'].
      ##    If this is the first call to this function for k' then this is exactly the initial GCOMP prediction.
      ##    Note that after we generate an update Qk_hat_star, we ****over-write**** current Qk_hat with its updated version
      Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]

      ## 4A. The model update. Univariate logistic regression (TMLE)
      if (Qk_idx == max_Qk_idx) {
        if (gvars$verbose) cat("Last targeting step for E(Y_d) with intercept only TMLE updates\n")
        ## Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]

        X <- as.matrix(qlogis(Qk_hat))
        newX <- as.matrix(qlogis(Qk_hat_all))
        colnames(X) <- colnames(newX) <- "offset"
        # TMLE.fit <- TMLE.updater.speedglm(Y = Qkplus1, X = X, newX = newX, obsWeights = wts)
        TMLE.fit <- TMLE.updater.glm(Y = Qkplus1, X = X, newX = newX, obsWeights = wts)
        Qk_hat_star_all <- TMLE.fit[["pred"]]
        Qk_hat_star <- predict(TMLE.fit[["fit"]], X)
        # TMLE.fit <- tmle.update(Qkplus1 = Qkplus1,
        #                         Qk_hat = Qk_hat,
        #                         IPWts = wts,
        #                         lower_bound_zero_Q = FALSE,
        #                         skip_update_zero_Q = self$skip_update_zero_Q)
        # TMLE_intercept <- TMLE.fit$TMLE_intercept
        # if (!is.na(TMLE_intercept) && !is.nan(TMLE_intercept)) {
        #   update.Qstar.coef <- TMLE_intercept
        # } else {
        #   update.Qstar.coef <- 0
        # }
        # Qk_hat_star_all2 <- plogis(qlogis(Qk_hat_all) + update.Qstar.coef)
        # max(Qk_hat_star_all2-Qk_hat_star_all)

        ## TO DO: REPLACE Qk_hat with Qk_hat_star (targeted version)
        EIC_i_t_calc <- wts * (Qkplus1 - Qk_hat_star)
        # EIC_i_t_calc <- wts * (Qkplus1 - Qk_hat)
        data$dat.sVar[use_subset_idx, ("EIC_i_t") := EIC_i_t_calc]

      ## 4B. The model update. Infinite dimensional epsilon (SDR)
      } else {
        ## Fitting SDR update with split-spec SL: Qkplus1 ~ offset(qlogis(Qk_hat)) + H[k-1] and weights 'wts'
        ## Grabbing covariates for H(k-1), where k-1 is determined by the row index offset Hk_row_offset
        ## Predict for all newly censored and new non-follows at k'
        ##    GENERALIZE THIS TO MAKE PREDICTIONS of E[Y_d|W'] for any new W'  BASED ON THIS MODEL UPDATE
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

        mfit <- iTMLE.updater.xgb(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts, params = self$reg$SDR_model)
        # mfit <- iTMLE.updater.glm(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # mfit <- TMLE.updater(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # mfit <- TMLE.updater.NULL(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        Qk_hat_star_all <- mfit[["pred"]]

        # ## GAM
        # # require('gam')
        # # form <- as.formula(paste0("Qkplus1 ~ offset(qlogis(Qk_hat))+", paste(names(obs_dat), collapse = "+")))
        # # mfit <- gam(
        # #   form,
        # #   data = cbind(Qkplus1 = Qkplus1, obs_dat),
        # #   y = Qkplus1,
        # #   family = binomial(),
        # #   # offset = offset(qlogis(Qk_hat)),
        # #   weights = wts)

        # # mfit <- gam.fit(
        # #   x = as.matrix(obs_dat),
        # #   y = Qkplus1,
        # #   smooth.frame = NA,
        # #   family = binomial(),
        # #   offset = qlogis(Qk_hat),
        # #   weights = wts)
        # # pred1.star = as.numeric(predict(mfit,type='response'))
        # ## GAM:
        # # Qk_hat_star_all <- as.numeric(predict(mfit, pred_dat[1:582,], type='link'))
        # # Qk_hat_star_all <- plogis(qlogis(Qk_hat_all) + Qk_hat_star_all)

      }

      # print("MSE of previous Qk_hat vs. upated Qk_hat: " %+% mean((Qk_hat_star_all-Qk_hat_all)^2))
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
        # if (gvars$verbose) cat("reached the last targeting iteration of the very last initial regression. Saving the final targeted prediction for E[Y_d]")
        private$probAeqa <- Qk_hat_star_all
      } else {
        data$dat.sVar[(self$subset_idx - 1), "Qkplus1" := Qk_hat_star_all]
      }
      invisible(self)
    }
  )
)