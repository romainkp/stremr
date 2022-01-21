## ---------------------------------------------------------------------
## R6 Class for Q-Learning (Iterative G-COMP or TMLE)
## ---------------------------------------------------------------------
##  R6 class for controlling the internal implementation of Q-learning functionality.
##  Supports sequential (recursive) G-computation and longitudinal TMLE.
##  Inherits from \code{ModelBinomial} R6 Class.
##
## @docType class
## @format An \code{\link{R6Class}} generator object
## @keywords R6 class
## @details
## \itemize{
##    \item{\code{regimen_names}} - \code{character}. Note used. For future pooling across regimens.
##    \item{\code{classify}} - ... \code{logical}
##    \item{\code{TMLE}} - ... \code{logical}
##    \item{\code{nIDs}} - ... \code{integer}
##    \item{\code{stratifyQ_by_rule}} - ... \code{logical}
##    \item{\code{lower_bound_zero_Q}} - ... \code{logical}
##    \item{\code{skip_update_zero_Q}} - ... \code{logical}
##    \item{\code{Qreg_counter}} - ... \code{integer}
##    \item{\code{t_period}} - ... \code{integer}
##    \item{\code{idx_used_to_fit_initQ}} - ... \code{integer}
## }
## @section Methods:
## \describe{
##   \item{\code{new(reg, ...)}}{...}
##   \item{\code{define.subset.idx(data, subset_vars, subset_exprs)}}{...}
##   \item{\code{define_idx_to_fit_initQ(data)}}{...}
##   \item{\code{define_idx_to_predictQ(data)}}{...}
##   \item{\code{fit(overwrite = FALSE, data, ...)}}{...}
##   \item{\code{Propagate_TMLE_fit(overwrite = TRUE, data, new.TMLE.fit, ...)}}{...}
##   \item{\code{predict(newdata, subset_idx, ...)}}{...}
##   \item{\code{predictStatic(data, g0, gstar, subset_idx)}}{...}
##   \item{\code{predictStochastic(data, g0, gstar, subset_idx, stoch_indicator)}}{...}
##   \item{\code{predictAeqa(newdata, ...)}}{...}
##   \item{\code{get.fits()}}{...}
## }
## @section Active Bindings:
## \describe{
##    \item{\code{wipe.alldat}}{...}
##    \item{\code{getfit}}{...}
##    \item{\code{getTMLEfit}}{...}
## }
## @export
ModelQlearn  <- R6Class(classname = "ModelQlearn",
  inherit = ModelBinomial,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    regimen_names = character(), ## for future pooling across regimens
    classify = FALSE,
    TMLE = TRUE,
    CVTMLE = FALSE,
    byfold_Q = FALSE,
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,   ## train only among those who are following the rule of interest?
    lower_bound_zero_Q = TRUE,
    skip_update_zero_Q = TRUE,
    Qreg_counter = integer(),    ## Counter for the current sequential Q-regression (min is at 1)
    Qreg_fail = FALSE,           ## Indicator that this Q model fit has failed (was unable to train and/or predict)
    all_Qregs_indx = integer(),
    t_period = integer(),
    keep_idx = FALSE,            ## keep current indices (do not remove them right after fitting)
    keep_model_fit = TRUE,       ## keep the model fit object for current Q_k
    idx_used_to_fit_initQ = NULL,
    fold_y_names = NULL,
    lwr = 0.0,                   ## lower bound for Q predictions
    upr = 1.0,                   ## upper bound for Q predictions
    maxpY = 1.0,                 ## max incidence P(Y=1|...) for rare-outcomes TMLE, only works with learners that can handle logistic-link with outcomes > 1
    TMLE_updater = NULL,         ## function to use for TMLE updates

    initialize = function(reg, ...) {
      super$initialize(reg, ...)

      self$reg <- reg
      self$outcome_type <- "quasibinomial"

      self$Qreg_counter <- reg$Qreg_counter
      self$all_Qregs_indx <- reg$all_Qregs_indx

      self$stratifyQ_by_rule <- reg$stratifyQ_by_rule
      self$lower_bound_zero_Q <- reg$lower_bound_zero_Q
      self$skip_update_zero_Q <- reg$skip_update_zero_Q
      self$t_period <- reg$t_period
      self$regimen_names <- reg$regimen_names
      self$maxpY <- reg$maxpY

      if (!is.character(reg$TMLE_updater)) stop("'TMLE_updater' must be a character name of an R function.")
      updater_fun <- get(reg$TMLE_updater)
      if (!is.function(updater_fun)) stop("'TMLE_updater' must resolve to an existing function name.")
      self$TMLE_updater <- updater_fun

      self$TMLE <- reg$TMLE
      self$CVTMLE <- reg$CVTMLE
      self$byfold_Q <- reg$byfold_Q

      self$keep_idx <- reg$keep_idx
      self$keep_model_fit <- reg$keep_model_fit

      if (gvars$verbose == 2) {print("initialized Q class"); reg$show()}

      invisible(self)
    },

    define.subset.idx = function(data, subset_vars, subset_exprs) { data$evalsubst(subset_vars, subset_exprs) },

    # Add all obs that were also censored at t and (possibly) those who just stopped following the rule
    define_idx_to_predictQ = function(data) { self$define.subset.idx(data, subset_exprs = self$subset_exprs) },

    define_idx_to_fit_initQ = function(data) {
      ## self$subset_vars is the name of the outcome var, checks for (!is.na(self$subset_vars))
      ## self$subset_exprs is the time variable for selecting certain rows
      subset_idx <- self$define.subset.idx(data, subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)

      ## excluded all censored observations:
      subset_idx <- intersect(subset_idx, which(data$uncensored))

      ## if stratifying by rule, exclude all obs who stopped following the rule at some point prior t or at t:
      if (self$stratifyQ_by_rule) subset_idx <- intersect(subset_idx, which(data$follow_rule))

      return(subset_idx)
    },

    ## **********************************************************************
    ## Iteration step for G-COMP
    ## (*) Update are performed only among obs who were involved in fitting the initial Q
    ## (*) Predicted outcome from the previous Seq-GCOMP iteration (t+1) was previously saved in current t row under "Qkplus1".
    ## (*) If person failed at t, he/she could not have participated in model fit at t+1.
    ## (*) Thus the outcome at row t for new failures ("Qkplus1") is always deterministically assigned to 1.
    ## **********************************************************************
    fit = function(overwrite = FALSE, data, Qlearn.fit, ...) { # Move overwrite to a field? ... self$overwrite

      prevQ_indx <- which(self$all_Qregs_indx %in% (self$Qreg_counter+1))
      if (length(prevQ_indx) > 0) {
        prevQ <- Qlearn.fit$getPsAsW.models()[[prevQ_indx]]
        prevQreg_fail <- prevQ$Qreg_fail
        if (self$byfold_Q) fold_y_names <- prevQ$fold_y_names
      } else {
        fold_y_names <- NULL
        prevQreg_fail <- FALSE
      }

      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      ## **********************************************************************
      ## TRAINING STEP TMLE/GCOMP: Fit the initial Q regression
      ## Select all obs at t who were uncensored & possibly were following the rule & had no missing outcomes (!is.na(self$subset_vars))
      ## Fit model using Q.kplus as the outcome to obtain the inital model fit for Q_k
      ## **********************************************************************
      self$subset_idx <- self$define_idx_to_fit_initQ(data)
      self$idx_used_to_fit_initQ <- self$subset_idx ## Save the subset used for fitting of the current initial Q[t] -> will be used for targeting
      nodes <- data$nodes
      self$n_obs_fit <- length(self$subset_idx)
      ## multiply the outcomes by delta-factor (for rare-outcomes adjustment with known upper bound = self$maxpY)
      data$rescaleNodes(subset_idx = self$subset_idx, nodes_to_rescale = self$outvar, delta = 1/self$maxpY)
      ## Fit model using (re-scaled) Q.kplus as the outcome to obtain the inital model fit for Q[t]:
      private$model.fit <- fit_single_regression(data = data, 
                                                 nodes = nodes, 
                                                 models = self$models, 
                                                 model_contrl = self$model_contrl, 
                                                 predvars = self$predvars, 
                                                 outvar = self$outvar, 
                                                 subset_idx = self$subset_idx,
                                                 outcome_type = self$outcome_type,
                                                 fold_y_names = fold_y_names)
      if (gvars$verbose) {
        print("Q.init model fit:"); try(print(private$model.fit))
      }
      ## Convert the outcome column back to its original scale (multpily by delta)
      data$rescaleNodes(subset_idx = self$subset_idx, nodes_to_rescale = self$outvar, delta = self$maxpY)
      self$is.fitted <- TRUE

      ## **********************************************************************
      ## FIRST PREDICTION STEP OF TMLE/GCOMP (Qk_hat = initial Q)
      ## Predict for subjects that were used for training initial Q
      ## E(Q_{k+1}|O(t)): Predict from conditional mean fit using the observed data (O) (no counterfactuals)
      ## TMLE: Qk_hat will be modified into Q^* with TMLE update
      ## GCOMP: Qk_hat serves no purpose and will be discarded
      ## **********************************************************************
      Qk_hat <- try(self$predict(data, subset_idx = self$idx_used_to_fit_initQ))
      # prediction in seq-Gcomp has failed
      if (inherits(Qk_hat, "try-error")) {
        warning("error during prediction step of GCOMP/TMLE, assigning NA to predicted Q values")
        print(Qk_hat)
        Qk_hat <- NA
        self$Qreg_fail <- TRUE
      } else {
        ## bound predictions
        Qk_hat[Qk_hat < self$lwr] <- self$lwr
        Qk_hat[Qk_hat > self$upr] <- self$upr
        ## rescale prediction from (0,1) range to its pre-specified bounded range (maximum probability constraint)
        Qk_hat <- Qk_hat*self$maxpY
      }

      ## TMLE update of initial prediction (will return Qk_hat if not doing TMLE updates)
      TMLE.fit <- self$TMLE_update(data, Qk_hat)

      ## **********************************************************************
      ## SECOND PREDICTION STEP OF TMLE/GCOMP
      ## 1. Predict for subjects that were used for training initial Q
      ## 2. Predict for subjects who were newly censored and just stopped following the rule
      ## E(Q_{k+1}|O^*(t)): Predict from conditional mean fit using the counterfactual (intervened) data (O^*)
      ## Set the intervention nodes (A(t),N(t)) to counterfactual values (A^*(t),N^*(t)), then predict
      ## TMLE: Qk_hat will be modified into Q^* with TMLE update
      ## GCOMP: Qk_hat serves no purpose and will be discarded
      ## **********************************************************************
      ## **********************************************************************
      ## ALGORITHM OUTLINE
      ## 1. Reset current exposures at t to their counterfactual values (by renaming and swapping the columns in input data)
      ##    note: this step is unnecessary when doing fitting only among rule-followers
      ## 3. Add observations that were censored at current t (+possibly stopped following the rule)
      ## 2. Predict for all subjects, under counterfactual regimen A(t)=A^*(t)
      ## 4. For stochastic g^*_t: Perform integration over the support of the stochastic intervention (as a weighted sum over g^*(a)=P(A^*(t)=a(t)))
      ## TO DO: add option for performing MC sampling to perform the same integration for stochastic g^*
      ## **********************************************************************
      interventionNodes.g0 <- data$interventionNodes.g0
      interventionNodes.gstar <- data$interventionNodes.gstar
      ## Determine which nodes are actually stochastic and need to be summed out:
      stoch_indicator <- data$define.stoch.nodes(interventionNodes.gstar)
      any_stoch <- sum(stoch_indicator) > 0
      ## For prediction need to add all obs that were also censored at t and (possibly) those who just stopped following the rule
      self$subset_idx <- self$define_idx_to_predictQ(data)

      if (!any_stoch) {
        pred_Qk <- self$predictStatic(data,
                                    g0 = interventionNodes.g0,
                                    gstar = interventionNodes.gstar,
                                    subset_idx = self$subset_idx,
                                    TMLE.fit = TMLE.fit)
      } else {
        # For all stochastic nodes, need to integrate out w.r.t. the support of each node
        pred_Qk <- self$predictStochastic(data,
                                        g0 = interventionNodes.g0,
                                        gstar = interventionNodes.gstar,
                                        subset_idx = self$subset_idx,
                                        stoch_indicator = stoch_indicator,
                                        TMLE.fit = TMLE.fit)
      }

      pred_Qk[pred_Qk < self$lwr] <- self$lwr
      pred_Qk[pred_Qk > self$upr] <- self$upr

      ## rescale the prediction in (0,1) range to its pre-specified bounded range (maximum probability constraint)
      pred_Qk <- pred_Qk*self$maxpY

      ## pred_Qk is the prediction of the target parameter (\psi_hat) at the current time-point k (where we already set A(k) to A^*(k))
      ## This is either the initial G-COMP Q or the TMLE targeted version of the initial Q
      ## This prediction includes all newly censored observations.
      ## When stratifying Q fits, these predictions will also include all observations who just stopped following the treatment rule at time point k.
      data$dat.sVar[self$subset_idx, ("Qk_hat") := pred_Qk]

      ## Set the outcome for the next Q-regression: put Q[t] in (t-1).
      ## only set the Qkplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
        newsubset_idx <- (self$subset_idx - 1)
        ## The outcomes of these IDs (row numbers) need to be replaced with fold-specific outcomes for each training set
        ## Next training fold model can check if some of its indices overlap with 'train_idx_to_replace'
        ## How do we pass these predictions and store them?
        ## train_idx_to_replace <- (self$idx_used_to_fit_initQ - 1)
        data$dat.sVar[newsubset_idx, ("Qkplus1") := pred_Qk]

        if (self$byfold_Q) {
          folds_seq <- seq_along(private$probA1_byfold)
          self$fold_y_names <- paste0("Qkplus1_f", folds_seq)
          ## initialize the values of fold-specific outcomes for next Q iteration, this is done for the very first Q regression:
          if (is.null(fold_y_names))
            data$dat.sVar[, (self$fold_y_names) := as.numeric(get(data$nodes$Ynode))]
          ## define split-specific / fold-specific Q predictions in the dataset as the outcomes for next Q regression  (separate column for each fold)
          for (split in folds_seq) {
            data$dat.sVar[newsubset_idx, paste0("Qkplus1_f", split) := private$probA1_byfold[[split]]]
          }
        }

        private$probA1 <- NULL

      } else {
        ## save prediction P(Q.kplus=1):
        private$probAeqa <- pred_Qk
        private$probA1 <- NULL
      }
      
      self$Qreg_fail <- prevQreg_fail | self$Qreg_fail
      
      # **********************************************************************
      # Wipe out all internal data -- DOES NOT REMOVE THE Q MODEL FIT OBJECTS
      self$wipe.alldat
      # **********************************************************************
      ## If we are planning on running iterative TMLE we will need the subsets used for fitting and predicting this Q
      if (!self$keep_idx) self$wipe.all.indices
      # **********************************************************************
      ## For CV, remove all grid & fold-specific model fits from gridisl
      ## Will still keep the best model fit that was re-trained on all data method=="cv"/"holdout"
      ## Don't need to store fold-specific models, since we already made the predictions.
      ## However, we may need them if we want to report the model fit stats
      if (!(self$model_contrl[["fit_method"]] %in% "none")) private$model.fit$wipe.allmodels
      # **********************************************************************

      invisible(self)
    },

    TMLE_update = function(data, init_Qk) {
      if (!self$TMLE) {

        return(NULL)

      } else if (self$TMLE) {

        ## The outcome (not re-scaled yet by 1/self$maxpY) that was used for fitting the initial Q at current time-point k:
        Qkplus1 <- data$dat.sVar[self$idx_used_to_fit_initQ, self$outvar, with = FALSE][[1]]
        ## re-scaled initial predictions (from initial outcome model)
        init_Qk <- init_Qk / self$maxpY
        ## re-scaled outcomes
        Y <- Qkplus1 / self$maxpY

        ## TMLE offset (log(x/[1-x])) is derived from the initial prediction of Q among ROWS THAT WERE USED TO FIT Q
        ## Thus, need to find which elements in predicted Q vector (init_Qk) where actually used for fitting the init Q
        # idx_for_fits_among_preds <- which(self$subset_idx %in% self$idx_used_to_fit_initQ)
        # init_Qk_fitted <- init_Qk[idx_for_fits_among_preds]
        ## Cumulative IPWeights for current t=k (from t=0 to k):
        wts_TMLE <- data$IPwts_by_regimen[self$idx_used_to_fit_initQ, "cum.IPAW", with = FALSE][[1]]

        ## TMLE update based on the IPWeighted logistic regression model with offset and intercept only:
        # newX <- X <- data.frame(intercept = rep(1L, length(init_Qk_fitted)), offset = qlogis(init_Qk_fitted))
        newX <- X <- data.frame(intercept = rep(1L, length(init_Qk)), offset = qlogis(init_Qk))

        ## todo: in case of failure call TMLE.updater.NULL, allow passing custom TMLE updaters
        xgb.params <- list(nrounds = 100, silent = 1, objective = "reg:logistic", booster = "gblinear", eta = 1, nthread = 1)

        TMLE.fit <- self$TMLE_updater(Y, X, newX, obsWeights = wts_TMLE, params = xgb.params)
        # TMLE.fit <- iTMLE.updater.xgb(Y, X, newX, obsWeights = wts_TMLE, params = xgb.params)
        # TMLE.fit <- TMLE.updater.speedglm(Y, X, newX, obsWeights = wts_TMLE)
        # TMLE.fit <- TMLE.updater.glm(Y, X, newX, obsWeights = wts_TMLE)
        # TMLE.fit <- TMLE.updater.NULL(Y, X, newX, obsWeights = wts_TMLE)

        ## Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        Qk_hat <- TMLE.fit$pred*self$maxpY
        ## evaluate the contribution to the EIC (TMLE variance estimation)
        # Qk_hat <- pred_Qk[idx_for_fits_among_preds]
        EIC_i_t_calc <- wts_TMLE * (Qkplus1 - Qk_hat)
        data$dat.sVar[self$idx_used_to_fit_initQ, ("EIC_i_t") := EIC_i_t_calc]

        return(TMLE.fit$fit)
      }
    },

    # **********************************************************************
    # Take a new TMLE fit and propagate it by first updating the Q(t) model and then updating the Q-model-based predictions for previous time-point
    # **********************************************************************
    Propagate_TMLE_fit = function(overwrite = TRUE, data, new.TMLE.fit, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitly asked
      # TMLE_intercept <- new.TMLE.fit$TMLE_intercept
      # if (!is.na(TMLE_intercept) && !is.nan(TMLE_intercept)) {
      #   update.Qstar.coef <- TMLE_intercept
      # } else {
      #   update.Qstar.coef <- 0
      # }

      ## Update the model predictions (Qk_hat) for initial Q[k] from GCOMP at time-point k.
      ## Based on TMLE update, the predictions now include ALL obs that are newly censored and just stopped following the rule at k:
      Qk_hat <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]

      if (!is.null(new.TMLE.fit)) {
        save <- private$probA1
        newdata <- data.frame(intercept = rep(1L, length(Qk_hat)))
        Qk_hat_updated <- predict(new.TMLE.fit, newdata = newdata, offset = qlogis(Qk_hat))
        # private$probA1 <- plogis(qlogis(private$probA1) + update.Qstar.coef)
        # head(cbind(private$probA1, save))
      } else {
        Qk_hat_updated <- Qk_hat
      }

      # Qk_hat_updated <- plogis(qlogis(Qk_hat) + update.Qstar.coef)

      # Save all predicted vals as Qk_hat[k] in row k or first target and then save targeted values:
      data$dat.sVar[self$subset_idx, "Qk_hat" := Qk_hat_updated]

      # Set the outcome for the next Q-regression: put Q[t] in (t-1), this will be overwritten with next prediction
      # only set the Qkplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
        data$dat.sVar[(self$subset_idx-1), "Qkplus1" := Qk_hat_updated]
        private$probA1 <- NULL
      } else {
        # save prediction P(Q.kplus=1) only for a subset in self$getsubset that was used for prediction
        private$probAeqa <- Qk_hat_updated
        private$probA1 <- NULL
      }

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # self$wipe.alldat
      # **********************************************************************
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw); Uses private$model.fit to generate predictions for data:
    predict = function(newdata, subset_idx, ...) {
      probA1 <- NULL
      probA1_byfold <- NULL
      ## For CV-TMLE need to use holdout predictions
      ## These are also referred to as predictions from all validation splits, or holdouts or out-of-sample predictions
      # holdout <- self$CVTMLE
      assert_that(self$is.fitted)
      model.fit <- private$model.fit

      if (missing(newdata) && !is.null(private$probA1)) {
        ## probA1 will be a one column data.table, hence we extract and return the actual vector of predictions:
        return(private$probA1)
      } else if (missing(newdata) && is.null(private$probA1)) {
        if (is(model.fit, "PredictionStack")) {
          probA1 <- gridisl::predict_SL(modelfit = model.fit, add_subject_data = FALSE, subset_idx = subset_idx, holdout = self$CVTMLE, verbose = gvars$verbose)
        } else if (is(model.fit, "Lrnr_base")) {
          probA1 <- model.fit$predict()
        } else {
          stop("model fit object is of unrecognized class (private$model.fit)")
        }

        ## probA1 will be a one column data.table, hence we extract and return the actual vector of predictions:
        private$probA1 <- probA1[[1]]
        return(private$probA1)

      } else {

        self$n <- newdata$nobs
        if (missing(subset_idx)) {
          subset_idx <- self$define.subset.idx(newdata, subset_exprs = self$subset_exprs)
        }

        if (!self$CVTMLE) {
          if (is(model.fit, "PredictionStack")) {
            probA1 <- gridisl::predict_SL(modelfit = model.fit,
                                          newdata = newdata,
                                          add_subject_data = FALSE,
                                          subset_idx = subset_idx,
                                          # use_best_retrained_model = TRUE,
                                          holdout = FALSE,
                                          verbose = gvars$verbose)
            if (self$byfold_Q) {
              probA1_byfold <- gridisl::predict_SL(modelfit = model.fit,
                                         newdata = newdata,
                                         add_subject_data = FALSE,
                                         subset_idx = subset_idx,
                                         byfold = TRUE,
                                         # use_best_retrained_model = TRUE,
                                         verbose = gvars$verbose)
              private$probA1_byfold <- probA1_byfold[["pred"]]
            }
          } else if (is(model.fit, "Lrnr_base")) {
            ## todo: allow passing the DataStorage object directly to task, seamlessly
            new_task <- sl3::sl3_Task$new(newdata$dat.sVar[self$subset_idx, ], covariates = self$predvars, outcome = self$outvar)
            probA1 <- model.fit$predict(new_task)
          } else {
            stop("model fit object is of unrecognized class (private$model.fit)")
          }

        } else {
          probA1 <- gridisl::predict_SL(modelfit = model.fit,
                                        newdata = newdata,
                                        add_subject_data = FALSE,
                                        subset_idx = self$idx_used_to_fit_initQ,
                                        # use_best_retrained_model = TRUE,
                                        holdout = TRUE,
                                        verbose = gvars$verbose)
          probA1[, ("idx") := self$idx_used_to_fit_initQ]

          newObs_idx <- setdiff(subset_idx, self$idx_used_to_fit_initQ)
          if (length(newObs_idx) > 0) {
            probA1_newObs <- gridisl::predict_SL(modelfit = model.fit,
                                                 newdata = newdata,
                                                 add_subject_data = FALSE,
                                                 subset_idx = newObs_idx,
                                                 # use_best_retrained_model = TRUE,
                                                 holdout = FALSE,
                                                 verbose = gvars$verbose)
            probA1_newObs[, ("idx") := newObs_idx]
            probA1 <- data.table::rbindlist(list(probA1, probA1_newObs))
          }

          setkeyv(probA1, "idx")
          probA1[, ("idx") := NULL]
        }

        ## if probA1 will be a one column data.table, we extract and return the actual vector of predictions:
        if (is.list(probA1) || is.data.table(probA1) || is.data.frame(probA1)) {
          probA1 <- probA1[[1]]
        }

        ## if one col matrix, extract the first column that is assumed to contain predictions
        if (is.matrix(probA1)) {
          probA1 <- probA1[,1]
        }

        ## check that predictions P(A=1 | dmat) exist for all obs
        if (any(is.na(probA1) & !is.nan(probA1))) {
          stop("some of the predicted probabilities during seq g-comp resulted in NAs, which indicates an error of a prediction routine")
        }
        private$probA1 <- probA1
        return(private$probA1)
      }
    },

    predictStatic = function(data, g0, gstar, subset_idx, TMLE.fit = NULL) {
      # ------------------------------------------------------------------------------------------------------------------------
      # Set current A's and N's to the counterfactual exposures in the data (for predicting Q):
      # ------------------------------------------------------------------------------------------------------------------------
      data$swapNodes(current = g0, target = gstar)

      # ------------------------------------------------------------------------------------------------------------------------
      # Predict based on counterfactual exposure settings
      # ------------------------------------------------------------------------------------------------------------------------
      gcomp.pred.res <- try(self$predict(data, subset_idx = subset_idx))

      # ------------------------------------------------------------------------------------------------------------------------
      # Reset back the observed exposure to A[t] (swap back by renaming columns)
      # ------------------------------------------------------------------------------------------------------------------------
      data$swapNodes(current = gstar, target = g0)

      # prediction in seq-Gcomp has failed
      if (inherits(gcomp.pred.res, "try-error")) {
        private$probA1 <- NA
        self$Qreg_fail <- TRUE
        warning("error during prediction step of GCOMP/TMLE, assigning NA to predicted Q values")
        print(gcomp.pred.res)
      }

      # ------------------------------------------------------------------------------------------------------------------------
      # Perform qlogit-linear TMLE update for predicted Q, if available
      # ------------------------------------------------------------------------------------------------------------------------
      if (!is.null(TMLE.fit)) {
        save <- private$probA1
        newdata <- data.frame(intercept = rep(1L, length(private$probA1)))
        private$probA1 <- predict(TMLE.fit, newdata = newdata, offset = qlogis(private$probA1))
        # private$probA1 <- plogis(qlogis(private$probA1) + update.Qstar.coef)
        # head(cbind(private$probA1, save))
      }

      invisible(return(private$probA1))
    },

    predictStochastic = function(data, g0, gstar, subset_idx, stoch_indicator, TMLE.fit = NULL) {
      # Dimensionality across all stochastic nodes:
      stoch_nodes_idx <- which(stoch_indicator)
      stoch_nodes_names <- names(stoch_indicator[stoch_nodes_idx])
      bit_list <- rep.int(list(c(0,1)), length(stoch_nodes_names))
      # Create a grid matrix, a single loop over the support of all nodes is a loop over the rows on the matrix
      all_vals_mat <- do.call(expand.grid, bit_list)
      d_all <- nrow(all_vals_mat)
      colnames(all_vals_mat) <- stoch_nodes_names
      # Save the probability of each stochastic intervention node
      stoch.probs <- data$dat.sVar[subset_idx, stoch_nodes_names, with = FALSE]
      colnames(stoch.probs) <- stoch_nodes_names

      # Loop over the grid mat; don't need to evaluate for everyone, just those obs that were used in prediction:
      stoch.probA1 <- 0
      for (i in 1:nrow(all_vals_mat)) {
        # Assign the values in all_vals_mat[i,] to stochastic nodes in newdata
        # modify data to assign a single value from the support of each stochastic node TO all observations
        # WARNING: THIS STEP IS IRREVERSIBLE, ERASES ALL CURRENT VALUES IN interventionNodes.g0[stoch_nodes_idx]:
        data$dat.sVar[subset_idx, (stoch_nodes_names) := all_vals_mat[i,]]

        # Predict using newdata, obtain probAeqa_stoch
        # Predict Prob(Q.init = 1) for all observations in subset_idx (note: probAeqa is never used, only private$probA1)
        # Obtain the initial Q prediction P(Q=1|...) for EVERYBODY (including those who just got censored and those who just stopped following the rule)
        probA1 <- self$predictStatic(data, g0 = g0, gstar = gstar, subset_idx = subset_idx, TMLE.fit = TMLE.fit)
        # Evaluate the joint probability vector for all_vals_mat[i,] for n observations (cumulative product) based on the probabilities from original column
        jointProb <- rep.int(1L, nrow(stoch.probs))
        for (stoch.node.nm in stoch_nodes_names) {
          IndNodeVal <- all_vals_mat[i, stoch.node.nm]
          stoch.prob <- stoch.probs[[stoch.node.nm]]
          stoch.prob <- (stoch.prob)^IndNodeVal * (1L-stoch.prob)^(1-IndNodeVal)
          jointProb <- jointProb * stoch.prob
          # put the probabilities back into input data:
          data$dat.sVar[subset_idx, (stoch.node.nm) := stoch.probs[[stoch.node.nm]]]
        }
        # Weight the current prediction by its probability
        probA1 <- probA1 * jointProb
        # probA1[subset_idx] <- probA1[subset_idx] * jointProb
        # Sum and keep looping
        stoch.probA1 <- stoch.probA1 + probA1
        # stoch.probA1 <- stoch.probA1 + probA1[subset_idx]
      }
      private$probA1 <- stoch.probA1
      # stoch.probA1 <- stoch.probA1
      # private$probA1[subset_idx] <- stoch.probA1
      return(invisible(private$probA1))
    },

    # Return the presaved prediction P(Q.kplus=1) only for a subset based on self$getsubset and private$probA1 (which has all n predictions)
    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      } else {
        stop("ModelQlearn$predictAeqa has nothing to return, previous predictions have been deleted")
      }
    },

    # Returns the object that contains the actual model fits (itself)
    get.fits = function() {
      model.fit <- self$getfit
      tmle.fit <- self$getTMLEfit
      return(list(model.fit, tmle.fit))
    }
  ),

  active = list(
    wipe.alldat = function() {
      private$probA1 <- NULL
      private$probA1_byfold <- NULL
      ## If we are planning on running iterative TMLE we will need the indicies used for fitting and predicting this Q
      if (!self$keep_idx) {
        self$wipe.all.indices
      }
      if (!self$keep_model_fit) {
        self$wipe.model.fit
      }
      return(invisible(self))
    },
    wipe.all.indices = function() {
      self$idx_used_to_fit_initQ <- NULL
      self$subset_idx <- NULL
      return(invisible(self))
    },
    wipe.model.fit = function() {
      private$model.fit <- NULL
      return(invisible(self))
    },
    wipe.probs = function() {
      private$probA1 <- NULL
      private$probAeqa <- NULL
      private$probA1_byfold <- NULL
      return(invisible(self))
    },
    getprobA1_byfold = function() { private$probA1_byfold },
    getfit = function() { private$model.fit },
    getTMLEfit = function() { private$TMLE.fit }
  ),

  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    TMLE.fit = list(NA),
    .outvar = NULL,
    probA1 = NULL,    # Predicted probA^s=1 conditional on Xmat
    probAeqa = NULL,   # Likelihood of observing a particular value A^s=a^s conditional on Xmat
    probA1_byfold = NULL
  )
)