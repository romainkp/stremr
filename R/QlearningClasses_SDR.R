

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
        # browser()
        ## Second loop for targeting:
        for (i in (1:k_i)) {
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
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,
    Qreg_counter = integer(), # Counter for the current sequential Q-regression (min is at 1)
    t_period = integer(),

    ## **********************************************************************
    ## NEW SDR TARGETING PROCEDURE
    ## **********************************************************************
    ## TO DO:
    ## **********************************************************************
    ## * Need to implement correct evaluation of product weights (Qk_idx)
    ##
    ## * When k=-1 (after the last run of the epsilon targeting regressions).
    ##    Last loop can be currently identified as (Qk_idx == max_Qk_idx). However, this targets the functional Q_K(W).
    ##    Hence, we still need to do ***ANOTHER LOOP*** that will target E(Q_K(W)).
    ##    This final does intercept-only logistic regression updates (essentially the usual TMLE?)
    ##    Need to figure out which values of i, Qk_idx, t_period correspond with last seqG-COMP regression
    ##    NEED TO FIGURE OUT HOW TO CALL THIS FINAL LOOP AND WHERE TO STORE THE PREDICTIONS.
    ##
    ## * Still need to figure out how to pass covariates from the current G-COMP fit of Q_k, the covariates need to be from k-1.
    ## * Add newly censored obs and newly non-followers to targeted predictions (see below)
    ## **********************************************************************

    ## This update is for current time-point k' aimed at targeting initial Q_k fit.
    ## Outcome: The same outcomes that were used to fit this initial Q_{k'}
    ## Predictors: All predictors used for fitting Q_{k-1}
    ## Weights: reverse.cum.IPAW - the product of g's from t=k to current k' (cumulative product in reverse time ordering?)
    ## Offset: qlogis(Q.kplus1) - current (initial) Q prediction (at k')
    fit_epsilon_Q_k = function(data, kprime_idx, Qk_idx, max_Qk_idx, ...) {
      if (self$all_Qregs_indx[kprime_idx] != self$Qreg_counter) stop("something terrible has happened")
      cat("...SDR targeting...\n")
      cat("Total number of time-points: " %+% max_Qk_idx, "\n")

      cat("Current loop to update Q_k index = " %+% self$all_Qregs_indx[Qk_idx], "\n")

      cat("Targeting Q_k' (kprime) at = " %+% self$Qreg_counter, "; time = " %+% self$t_period, "\n")

      ## above is the same as self$all_Qregs_indx[kprime_idx]
      cat("Current k' (kprime) = " %+% kprime_idx, "\n")

      browser()

      # TMLE.intercept <- new.TMLE.fit$TMLE.intercept
      # if (!is.na(TMLE.intercept) && !is.nan(TMLE.intercept)) {
      #   update.Qstar.coef <- TMLE.intercept
      # } else {
      #   update.Qstar.coef <- 0
      # }
      self$t_period
      self$Qreg_counter

      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs

      ## use only the observations that participated in fitting of the initial Q_{k'} (current time-point is k')
      use_subset_idx <- self$idx_used_to_fit_initQ

      print("length(use_subset_idx): " %+% length(use_subset_idx))
      print("length(self$subset_idx): " %+% length(self$subset_idx))

      ## Weights: will define new column reverse.cum.IPAW
      wts <- OData$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]
      # wts <- data$IPwts_by_regimen[use_subset_idx, "reverse.cum.IPAW", with = FALSE][[1]]

      ## **TARGETTED** prediction from previous step at k'+1 (targeted towards Q_k).
      ## This is used as the outcome for the current regression update.
      ## INITIALIZED TO BE ALL Y AT FIRST
      Qstarkprime <- data$dat.sVar[use_subset_idx, "Qstarkprime", with = FALSE][[1]]

      ## Our offset is the Q that was previously targeted towards estimation of Q_k, now re-targeted towards Q_k'
      initQ <- data$dat.sVar[use_subset_idx, "Q.kplus1", with = FALSE][[1]]

      ## covariates \bar{H}(k-1)
      ## somehow we need to select from the same set of covariates only those people who are currently at risk

      require('xgboost')
      param <- list("objective" = "reg:logistic", "booster" = "gbtree", "nthread" = 4)
      obs_dat <- data$dat.sVar[use_subset_idx, self$predvars, with = FALSE]
      obs_dat[, CVD := as.numeric(CVD)]
      xgb_dat <- xgb.DMatrix(as.matrix(obs_dat), label = Qstarkprime)
      setinfo(xgb_dat, "base_margin", qlogis(initQ))
      setinfo(xgb_dat, "weight", wts)
      mfit <- xgb.train(params=param, data=xgb_dat, nrounds=10)

      ## NEED TO OBTAIN PREDICTIONS IN SAME WAY AS IN QlearnModel$fit, i.e.,
      ## ADD ALL NEWLY CENSORED OBSERVATIONS AND NEW NON-FOLLOWERS.
      ## THE wts used no longer play any role, but offsets need to be re-evaluated for new subset
      ## The predictors used are still all predictors from time point k-1
      ## (+ add all predictors from k-1 for all newly-censored non-followers for current time k').
      pred1.star = predict(mfit, xgb_dat)

      # require('gam')
      # mfit <- gam.fit(
      #   x=as.matrix(data$dat.sVar[use_subset_idx,self$predvars]),
      #   y = Qstarkprime,
      #   family=binomial(),
      #   offset = offset(qlogis(initQ)),
      #   weights=wts)
      # pred1.star = as.numeric(predict(mfit,type='response'))

      # Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
      # Q.kplus1 <- data$dat.sVar[self$subset_idx, "Q.kplus1", with = FALSE][[1]]
      # Q.kplus1.new <- plogis(qlogis(Q.kplus1) + update.Qstar.coef)

      # Save all predicted vals as Q.kplus1[t] in row t or first target and then save targeted values:
      # data$dat.sVar[self$subset_idx, "Q.kprime.2nd" := Q.kplus1.new]

      # Set the outcome for the next Q-regression: put Q[t] in (t-1), this will be overwritten with next prediction
      # only set the Q.kplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (kprime_idx < Qk_idx) {
        data$dat.sVar[(self$subset_idx - 1), "Qstarkprime" := pred1.star]
      } else {
        ## kprime_idx == Qk_idx means we reached the actual current Q_k that we were targeting,
        ## hence save it in the column where we store our global and final Q^* predictions
        data$dat.sVar[self$subset_idx, "Q.kplus1" := pred1.star]
      }
      invisible(self)
    },

    # Returns the object that contains the actual model fits (itself)
    get.fits = function() {
      model.fit <- self$getfit
      tmle.fit <- self$getTMLEfit
      return(list(model.fit, tmle.fit))
    }
  )
)