
## Split-specific predictions for validation sets
## Uses offsets that are also allowed to be split specific
cv_validation_preds <- function(fold, split_preds_Qk_hat) {
  ## These will be automatically defined in the calling frame of this function
  ## when the cross-validator that calls cv_split_preds()
  v <- origami::fold_index()
  train_idx <- origami::training()
  valid_idx <- origami::validation()
  v_Qk_hat <- as.numeric(split_preds_Qk_hat[["splitQk"]][[v]][valid_idx, ])
  list(v_Qk_hat = v_Qk_hat, valid_idx = valid_idx)
}

## Split-specific predictions for new data (for all nrow(data).
## Uses offsets that could might be also split-specific
cv_split_preds <- function(fold, data, fits_Qk, split_preds_Qk_hat, use_full = FALSE) {
  ## These will be automatically defined in the calling frame of this function
  ## when the cross-validator that calls cv_split_preds()
  v <- origami::fold_index()
  train_idx <- origami::training()
  valid_idx <- origami::validation()

  if ((all.equal(train_idx, valid_idx) == TRUE) || use_full) {
    # We're in final resubstitution call, so let's use full Q and g fits
    splitQk_fit <- fits_Qk$fullFit
  } else {
    # Split-specific Super Learners
    splitQk_fit <- fits_Qk$foldFits[[v]]
  }

  ## split-specific predictions for new data
  ## this may include new observations, e.g., extrapolating for newly censored
  ## The offsets could be split specific if there were previous SDR targeting steps
  new_data <- as.matrix(data)
  if (!is.null(split_preds_Qk_hat)) {
    new_data[, "offset"] <- qlogis(as.numeric(split_preds_Qk_hat[["splitQk"]][[v]]))
  }
  splitQk <- predict(splitQk_fit, newdata = new_data)[["pred"]]
  list(splitQk = splitQk)
}

## Ignores Y and instead uses Z generated from a split-specific
## training set generated using cv_split_preds X are the nodes to base the rule on
split_cv_SL <- function(fold, Y, X,
                        SL.library,
                        family,
                        obsWeights,
                        id,
                        use_full = FALSE,
                        split_preds_Qk_hat,
                        subset_Qk_hat,
                        split_preds_Qkplus1,
                        subset_Qplus1_newQ,
                        subset_Y, ...) {

  v <- origami::fold_index()
  train_idx <- origami::training()
  valid_idx <- origami::validation()

  ## -------------------------------------------------------------------------------------------------------
  ## ****** THIS GIVES US ACCESS TO SPLIT SPECIFIC PREDICTIONS (N in TOTAL) FROM SL trained on FOLD v ****
  ## -------------------------------------------------------------------------------------------------------
  # print(names(split_preds_Qk_hat))
  # print(names(split_preds_Qkplus1))
  # length(split_preds_Qk_hat[["splitQk"]])

  ## NEED TO USE NON-EXTRAPOLATED split-specific offsets (which are actually being used now for targeting)
  ## THIS LIST RECORDS ALL SPLIT-SPECIFIC model predictions (including extrapolated)
  # length(split_preds_Qk_hat[["splitQk"]][[1]])
  ## ******* Need to find out what indices in self$subset_idx ARE SPLIT-SPECIFIC (NOT NEW)
  ## ******* Those are the only outcomes that need to be replaced *******
  # length(split_preds_Qkplus1[["splitQk"]][[1]])

  ## This replaces the usual Y passed down to split_cv_SL() with FOLD-SPECIFIC predictions of Y (Yhat).
  ## Each Yhat is a vector of N predictions from the best model (SL) that was trained only on train_idx (v fold).

  ## WHAT ARE THE obsWeights? Are these the SL weights for each algorithm???
  # obsWeights <- obsWeights * split_preds$K[[v]]

  ## ADD NEW FAILURES (NEW OBSERVATIONS THAT WERE NOT USED FOR PREDICTING / FITTING Qkplus1 UNTIL NOW)
  if (!is.null(split_preds_Qkplus1)) {
    Y[subset_Y] <- as.numeric(split_preds_Qkplus1[["splitQk"]][[v]][subset_Qplus1_newQ, ])
  }

  ## ADD NEW SS-OFFSETS (ONLY WHEN AVAILABLE)
  ## NOTE THAT THESE NEED TO BE CONVERTED TO LOGIT-LINEAR SCALE!!!!
  if (!is.null(split_preds_Qk_hat)) {
    X[, "offset"] <- qlogis(as.numeric(split_preds_Qk_hat[["splitQk"]][[v]][subset_Qk_hat, ]))
  }

  origami::cv_SL(fold, Y, X, SL.library, family, obsWeights, id, ...)
}

## ---------------------------------------------------------------------
## R6 Class for Split-Specific Sequentially Double Robustness Targeting Procedure
## Internal implementation of Q-learning functionality.
## Inherits from \code{ModelQlearn} R6 Class.
## ---------------------------------------------------------------------
SplitCVSDRModelQlearn  <- R6Class(classname = "SplitCVSDRModelQlearn",
  inherit = SDRModelQlearn,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    split_preds_Qk_hat = NULL,
    SDR_SL_fit = NULL,

    ## Split-Specific SL updater
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
        cat("Regression Qk_idx order = " %+% self$all_Qregs_indx[Qk_idx], "\n")
        cat("Using the covariate space at t index = " %+% (self$all_Qregs_indx[kprime_idx]+Hk_row_offset), "\n")
      }

      ## 1. Weights: defined new column of cumulative weights where cumulative product starts at t = Qk_idx (k), rather than t = 0:
      wts <- data$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      ## 2a. Outcome: **TARGETED** prediction of the previous step k'+1.
      Qkplus1 <- data$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]
      ## 3. Offset:
      Qk_hat <- data$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]

      ## SS-SL will use split-specific offsets (Qhat) from this run and split-specific Y (Qkplus1) from the previous Q init / update (k+1)
      ## ** is.null(QModel_Qkplus1) means that there are no split-specific Y's available (use the actual observed Qkplus1 / Y)
      ##    If is.null(QModel_Qkplus1) then ALL Qkplus1 are the observed ones (no model predictions were used for assigning Qkplus1),
      ##    so no subsetting / replacing of Y is necessary
      ## ** is.null(split_preds_Qk_hat) means that no SS-offsets are available (first targeting step, only init Q fitted before)

      ## ** Subsets of subset_idx that need to be utilized during fitting, these are used for SS-offsets
      ##    IDs in split_preds_Qk_hat that need to be used for replacing offset in X with SS offsets (Qk_hat) for SS-targeting
      ##    ******** TO DO: WILL NEED SS-OFFSETS FOR prediction as well (extrapolating to new obs) ********
      subset_Qk_hat <- which(self$subset_idx %in% use_subset_idx)

      ## ** IDs in Y that need to be used to replace the original Y (Qkplus1) with split specific outcomes
      ##    This requires handling new observations (for which Qkplus1 is not a model prediction, but a deterministic value)
      subset_Y <- NULL
      subset_Qplus1_newQ <- NULL

      ## Previous round of Q fitting (k+1 time-point) is saved in QModel_Qkplus1 arg
      if (!is.null(QModel_Qkplus1)) {
        ## THESE ARE IDs in QModel_Qkplus1 that ARE ACTUALLY BEING USED FOR MODELING CURRENT Q.plus1
        subset_Qplus1_newQ <- which((QModel_Qkplus1$subset_idx - 1) %in% use_subset_idx)
        # data$dat.sVar[(QModel_Qkplus1$subset_idx - 1)[subset_Qplus1_newQ], ]
        ## THESE ARE IDs in Y that need to be replaced (not failures, i.e., why is not deterministic), must MATCH TO ABOVE IN LENGTH
        subset_Y <- which(use_subset_idx %in% (QModel_Qkplus1$subset_idx-1))
        # data$dat.sVar[use_subset_idx[subset_Y], ]
      }

      ## 4A. The model update. Univariate logistic regression (TMLE)
      if (Qk_idx == max_Qk_idx) {
        if (gvars$verbose) cat("Last targeting step for E(Y_d) with intercept only TMLE updates\n")
        ## Updated model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]

        X <- as.matrix(qlogis(Qk_hat))
        newX <- as.matrix(qlogis(Qk_hat_all))
        colnames(X) <- colnames(newX) <- "offset"

        ## This will use validation (out-of-sample) Qkplus1 & Qk_hat (offset) for the TMLE updates
        TMLE.fit <- try(TMLE.updater.speedglm(Y = Qkplus1, X = X, newX = newX, obsWeights = wts))
        if (inherits(TMLE.fit, "try-error")) { # TMLE update w/ speedglm failed
          TMLE.fit <- try(TMLE.updater.glm(Y = Qkplus1, X = X, newX = newX, obsWeights = wts))
        }
        Qk_hat_star_all <- TMLE.fit[["pred"]]
        Qk_hat_star <- predict(TMLE.fit[["fit"]], X)
        print("TMLE intercept update: " %+% TMLE.fit[["fit"]][["object"]][["coef"]])

        ## This will USE validation (out-of-sample) Qkplus1 & Qk_hat for the EIC estimates:
        ## REPLACE Qk_hat with Qk_hat_star (targeted version)
        EIC_i_t_calc <- wts * (Qkplus1 - Qk_hat_star)
        # EIC_i_t_calc <- wts * (Qkplus1 - Qk_hat)
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
        folds <- data$make_origami_fold_from_column(use_subset_idx + Hk_row_offset)

        # data$dat.sVar[self$subset_idx + Hk_row_offset,]
        pred_dat <- data$dat.sVar[self$subset_idx + Hk_row_offset, epsilon_predvars, with = FALSE]
        Qk_hat_all <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
        pred_dat[, ("offset") := qlogis(Qk_hat_all)]
        folds_pred_dat <- data$make_origami_fold_from_column(self$subset_idx + Hk_row_offset)
        # sort(unlist(lapply(folds_pred_dat, "[[", "validation_set")))

        # SL.library <- c( "TMLE.updater.NULL")
        SL.library <- c("TMLE.updater.NULL",
                      "TMLE.updater.glm",
                      "iTMLE.updater.glm",
                      "iTMLE.updater.xgb",
                      "iTMLE.updater.xgb.delta1",
                      "iTMLE.updater.xgb.delta2",
                      "iTMLE.updater.xgb.delta3",
                      "iTMLE.updater.xgb.delta4"
                      )

        library("abind")
        # SDR_SL_fit <- origami::origami_SuperLearner(folds = folds,
        #                                                    Y = Qkplus1,
        #                                                    X = as.matrix(obs_dat),
        #                                                    family = quasibinomial(),
        #                                                    obsWeights = wts,
        #                                                    SL.library = SL.library,
        #                                                    params = self$reg$SDR_model)

        ## Fit the Split-Specific SuperLearner.
        ## The function split_cv_SL() does all the work, taking the input data (X,Y) and
        ## replacing the training set Y[i] and offset[i] with split-specific outcomes Y^s[i] and offset^s[i]
        ## The rest of it goes on just like in the usual SuperLearner.
        ## The split-specific Y^s are obtained from split_preds_Qkplus1
        ##  (predictions Qhat from SS-SL for time-point t+1)
        ## The split-specific offset^s are obtained from split_preds_Qk_hat
        ##  (Qhat from the fit of the initial SL for time-point t or the targeted split-specific Qhat)
        SDR_SL_fit <- origami::origami_SuperLearner(folds = folds,
                                                    Y = Qkplus1,
                                                    X = as.matrix(obs_dat),
                                                    split_preds_Qk_hat = self$split_preds_Qk_hat,
                                                    subset_Qk_hat = subset_Qk_hat,
                                                    split_preds_Qkplus1 = QModel_Qkplus1$split_preds_Qk_hat,
                                                    subset_Qplus1_newQ = subset_Qplus1_newQ,
                                                    subset_Y = subset_Y,
                                                    family = quasibinomial(),
                                                    obsWeights = wts,
                                                    SL.library = SL.library,
                                                    params = self$reg$SDR_model,
                                                    cvfun = split_cv_SL)
        self$SDR_SL_fit <- SDR_SL_fit

        if (gvars$verbose) {
          print("Split CV SDR_SL_fit: "); print(SDR_SL_fit)
        }

        # names(SDR_SL_fit)
        # SDR_SL_fit_simple[["Z"]]
        # SDR_SL_fit[["Z"]]
        # SDR_SL_fit[["valY"]]
        # SDR_SL_fit[["valWeights"]]

        ## Slit-Specitic SuperLearner preds
        ## list of nrow(pred_dat) preds for each fold (including extrapolating preds for new obs)
        ## This also uses the offsets for newdata,
        ## if there were prior targeting steps, these will be split-specific offsets.
        ## All the split-specific predictions are being done here.
        ## This loops over each fold, then calls predict() on the fold-specific SL fits, passing newdata=data
        ## So each fold-specific SuperLearner fit does predictions for all observations in pred_dat.
        ## The results is a list of predictions: number of list items = number of folds.
        split_preds_Qk_hat <- origami::cross_validate(cv_split_preds,
                                                      folds,
                                                      data = pred_dat,
                                                      fits_Qk = SDR_SL_fit,
                                                      split_preds_Qk_hat = self$split_preds_Qk_hat,
                                                      .combine = FALSE)
        self$split_preds_Qk_hat <- split_preds_Qk_hat

        ## Validation (out-of-sample) predictions for the split-specific SuperLearner (n predictions by combining all validation folds in newdata)
        ## This function takes the list of previous split-spec predictions (above call to cross_validate) and
        ## extracts the validation set predictions from each fold. It then stacks these predictions together to obtain n validation predictions.
        Qk_hat_star_all <- origami::cross_validate(cv_validation_preds,
                                                   folds_pred_dat,
                                                   split_preds_Qk_hat = split_preds_Qk_hat,
                                                   .combine = TRUE)
        Qk_hat_star_all <- as.data.table(Qk_hat_star_all)
        setkeyv(Qk_hat_star_all, cols = "valid_idx")
        Qk_hat_star_all <- Qk_hat_star_all[["v_Qk_hat"]]

        ## Extract prediction for new data for regular SuperLearner (not split-specific, SL fit used ALL data)
        # Qk_hat_star_all <- as.numeric(predict(SDR_SL_fit, as.matrix(pred_dat))[["pred"]])

        # mfit <- iTMLE.updater.xgb(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts, params = self$reg$SDR_model)
        # mfit <- iTMLE.updater.glm(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # mfit <- TMLE.updater(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # mfit <- TMLE.updater.NULL(Y = Qkplus1, X = as.matrix(obs_dat), newX = as.matrix(pred_dat), obsWeights = wts)
        # cbind(Qk_hat_star_all, Qk_hat_star_all_2, Qk_hat_star_all - Qk_hat_star_all_2)
        # data_compare_fits <- data.frame(cbind(data$dat.sVar[self$subset_idx, "Qk_hat"], v_Qk_hat_star = Qk_hat_star_all, SL_Qk_hat_star = Qk_hat_star_all_2))
        # head(data_compare_fits)
      }

      # Over-write the old predictions with new model updates as Qk_hat[k'] in row [k']:
      data$dat.sVar[self$subset_idx, "Qk_hat" := Qk_hat_star_all]
      # Qk_hat_star_all-data$dat.sVar[self$subset_idx, Qk_hat]
      # print("MSE of previous Qk_hat vs. upated Qk_hat: " %+% mean((Qk_hat_star_all-Qk_hat_all)^2))

      ## ************** NOTE ******************
      ##  The shift below will be wrong when we get to the very last iteration of the very last update (LTMLE)
      ##  There will be nowhere to write a new update, since we will be at a row t==0.
      ##  This happens when: (Qk_idx == max_Qk_idx) && (Qk_idx == kprime_idx)
      ## **************************************
      ## Set the outcome for the next Q-regression: put the updated Q[k'] in (k'-1)

      ## (kprime_idx == Qk_idx) means that we have reached the final update for current Q[k] that we were targeting,
      ## Save this update in row [k'-1] = [k-1], that's where it will be picked up by next loop (or initial est step for Q[k-1]).
      # if ((Qk_idx == max_Qk_idx) && (Qk_idx == kprime_idx)) {
      #   if (gvars$verbose) cat("reached the last targeting iteration of the very last initial regression. Saving the final targeted prediction for E[Y_d]")
      #   private$probAeqa <- Qk_hat_star_all
      # } else {
      #   # data$dat.sVar[self$subset_idx, ]
      #   data$dat.sVar[(self$subset_idx - 1), "Qkplus1" := Qk_hat_star_all]
      # }
      # data$dat.sVar[(self$subset_idx - 1),]

      if ((Qk_idx == max_Qk_idx) && (self$Qreg_counter > 1)) {
        # if (gvars$verbose)
        # cat("reached the last targeting iteration of the very last initial regression. Saving the final targeted prediction for E[Y_d]")
        data$dat.sVar[(self$subset_idx - 1), "Qkplus1" := Qk_hat_star_all]
      } else {
        private$probAeqa <- Qk_hat_star_all
      }
      invisible(self)
    }
  )
)