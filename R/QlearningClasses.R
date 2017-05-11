tmle.update <- function(Qkplus1, Qk_hat, IPWts, lower_bound_zero_Q = TRUE, skip_update_zero_Q = TRUE) {
  QY.star <- NA
  if (sum(abs(IPWts)) < 10^-9) {
    update.Qstar.coef <- 0
    # if (gvars$verbose)
    message("GLM TMLE update cannot be performed since all IP-weights are exactly zero, setting epsilon = 0")
    warning("GLM TMLE update cannot be performed since all IP-weights are exactly zero, setting epsilon = 0")
  } else if ((sum(Qkplus1[IPWts > 0]) < 10^-5) && skip_update_zero_Q) {
    update.Qstar.coef <- 0
  } else {
    #************************************************
    # TMLE update via weighted univariate ML (espsilon is intercept)
    #************************************************
    if (lower_bound_zero_Q) {
      Qkplus1[Qkplus1 < 10^(-4)] <- 10^(-4)
      Qk_hat[Qk_hat < 10^(-4)] <- 10^(-4)
    }
    off_TMLE <- qlogis(Qk_hat)
    # off_TMLE[off_TMLE == Inf] <- 20
    # off_TMLE[off_TMLE == -Inf] <- -20
    off_TMLE[off_TMLE >= 20] <- 20
    off_TMLE[off_TMLE <= -11] <- -11.0

    m.Qstar <- try(speedglm::speedglm.wfit(X = matrix(1L, ncol = 1, nrow = length(Qkplus1)),
                                          y = Qkplus1, weights = IPWts, offset = off_TMLE,
                                          # method=c('eigen','Cholesky','qr'),
                                          family = quasibinomial(), trace = FALSE, maxit = 1000),
                  silent = TRUE)

    if (inherits(m.Qstar, "try-error")) { # TMLE update failed
      # if (gvars$verbose)
      message("GLM TMLE update has failed, setting epsilon = 0")
      warning("GLM TMLE update has failed, setting epsilon = 0")
      update.Qstar.coef <- 0
    } else {
      update.Qstar.coef <- m.Qstar$coef
    }
  }

  if (gvars$verbose == 2) cat("...TMLE epsilon update = ", update.Qstar.coef, "\n")

  fit <- list(TMLE_intercept = update.Qstar.coef)
  class(fit)[2] <- "tmlefit"
  if (gvars$verbose) print("tmle update: " %+% update.Qstar.coef)
  return(fit)
}

## ---------------------------------------------------------------------
## R6 Class for Q-Learning (Iterative G-COMP or TMLE)
## ---------------------------------------------------------------------
##  R6 class for controlling the internal implementation of Q-learning functionality.
##  Supports sequential (recursive) G-computation and longitudinal TMLE.
##  Inherits from \code{BinaryOutcomeModel} R6 Class.
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
QlearnModel  <- R6Class(classname = "QlearnModel",
  inherit = BinaryOutcomeModel,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    regimen_names = character(), # for future pooling across regimens
    classify = FALSE,
    TMLE = TRUE,
    CVTMLE = FALSE,
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,
    lower_bound_zero_Q = TRUE,
    skip_update_zero_Q = TRUE,
    Qreg_counter = integer(), # Counter for the current sequential Q-regression (min is at 1)
    all_Qregs_indx = integer(),
    t_period = integer(),
    keep_idx = FALSE,         # keep current indices (do not remove them right after fitting)
    idx_used_to_fit_initQ = NULL,

    initialize = function(reg, ...) {
      super$initialize(reg, ...)

      self$reg <- reg

      self$Qreg_counter <- reg$Qreg_counter
      self$all_Qregs_indx <- reg$all_Qregs_indx

      self$stratifyQ_by_rule <- reg$stratifyQ_by_rule
      self$lower_bound_zero_Q <- reg$lower_bound_zero_Q
      self$skip_update_zero_Q <- reg$skip_update_zero_Q
      self$t_period <- reg$t_period
      self$regimen_names <- reg$regimen_names
      self$TMLE <- reg$TMLE
      self$CVTMLE <- reg$CVTMLE
      self$keep_idx <- reg$keep_idx

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

    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      # **********************************************************************
      # FITTING STEP OF Q-LEARNING
      # Select all obs at t who were uncensored & possibly were following the rule & had no missing outcomes (!is.na(self$subset_vars))
      # **********************************************************************
      self$subset_idx <- self$define_idx_to_fit_initQ(data)
      # Save the subset used for fitting of the current initial Q[t] -> will be used for targeting
      self$idx_used_to_fit_initQ <- self$subset_idx

      nodes <- data$nodes
      self$n_obs_fit <- length(self$subset_idx)

      # Fit model using Q.kplus as the outcome to obtain the inital model fit for Q[t]:
      private$model.fit <- fit_single_regression(data, nodes, self$models, self$model_contrl, self$predvars, self$outvar, self$subset_idx)
      self$is.fitted <- TRUE

      # **********************************************************************
      # PREDICTION STEP OF Q-LEARNING
      # Q prediction for everyone (including those who just got censored and those who just stopped following the rule)
      # **********************************************************************
      interventionNodes.g0 <- data$interventionNodes.g0
      interventionNodes.gstar <- data$interventionNodes.gstar

      # Determine which nodes are actually stochastic and need to be summed out:
      stoch_indicator <- data$define.stoch.nodes(interventionNodes.gstar)
      any_stoch <- sum(stoch_indicator) > 0

      # For prediction need to add all obs that were also censored at t and (possibly) those who just stopped following the rule
      self$subset_idx <- self$define_idx_to_predictQ(data)

      ## **********************************************************************
      ## Below prediction step will:
      ## **********************************************************************
      ## 1. reset current exposures at t to their counterfactual values (e.g., by renaming the columns)
      ## note: this step is unnecessary when doing fitting only among rule-followers
      ## self$binomialModelObj$replaceCols(data, self$subset_idx, data$nodes$, , ...)
      ## 2. predict for all obs who were used for fitting t (re-using the same design mat used for fitting)
      ## probAeqa <- self$predictAeqa(...)
      ## self$predict(...)
      ## 3. add observations that were censored at current t to the subset and do predictions for them
      ## ------------------------------------------------------------------------------------------------------------------------
      ## 4. possibly intergrate out over the support of the stochastic intervention (weighted sum of P(A(t)=a(t)))
      ## ------------------------------------------------------------------------------------------------------------------------
      ## *** NEED TO ADD: do MC sampling to perform the same integration for stochastic g^*
      ## ------------------------------------------------------------------------------------------------------------------------
      if (!any_stoch) {
        probA1 <- self$predictStatic(data,
                                    g0 = interventionNodes.g0,
                                    gstar = interventionNodes.gstar,
                                    subset_idx = self$subset_idx)
      } else {
        # For all stochastic nodes, need to integrate out w.r.t. the support of each node
        probA1 <- self$predictStochastic(data,
                                        g0 = interventionNodes.g0,
                                        gstar = interventionNodes.gstar,
                                        subset_idx = self$subset_idx,
                                        stoch_indicator = stoch_indicator)
      }

      iQ_all <- probA1

      ## print("initial mean(Qkplus1) for ALL obs at t=" %+% self$t_period %+% ": " %+% round(mean(iQ_all), 4))

      ## **********************************************************************
      ## Iteration step for G-COMP
      ## (*) Update are performed only among obs who were involved in fitting the initial Q (self$idx_used_to_fit_initQ)
      ## (*) Predicted outcome from the previous Seq-GCOMP iteration (t+1) was previously saved in current t row under "Qkplus1".
      ## (*) If person failed at t, he/she could not have participated in model fit at t+1.
      ## (*) Thus the outcome at row t for new failures ("Qkplus1") is always deterministically assigned to 1.
      ## **********************************************************************

      ## The outcome that was used for fitting the initial Q at current time-point k:
      Qkplus1 <- data$dat.sVar[self$idx_used_to_fit_initQ, "Qkplus1", with = FALSE][[1]]
      # data$dat.sVar[self$idx_used_to_fit_initQ, "Qkplus1" := eval(as.name("Qkplus1"))]

      if (self$TMLE) {
        ## TMLE offset (log(x/[1-x])) is derived from the initial prediction of Q among ROWS THAT WERE USED TO FIT Q
        ## Thus, need to find which elements in predicted Q vector (probA1) where actually used for fitting the init Q
        idx_for_fits_among_preds <- which(self$subset_idx %in% self$idx_used_to_fit_initQ)
        Qk_hat <- probA1[idx_for_fits_among_preds]
        # off_TMLE <- qlogis(Qk_hat)
        # print("initial mean(Qkplus1) among fitted obs only at t=" %+% self$t_period %+% ": " %+% round(mean(iQ_all), 4))

        ## Cumulative IPWeights for current t=k (from t=0 to k):
        wts_TMLE <- data$IPwts_by_regimen[self$idx_used_to_fit_initQ, "cum.IPAW", with = FALSE][[1]]

        ## TMLE update based on the IPWeighted logistic regression model with offset and intercept only:
        ## tmle update will error out if some predictions are exactly 0
        private$TMLE.fit <- tmle.update(Qkplus1 = Qkplus1,
                                        Qk_hat = Qk_hat,
                                        IPWts = wts_TMLE,
                                        lower_bound_zero_Q = self$lower_bound_zero_Q,
                                        skip_update_zero_Q = self$skip_update_zero_Q)

        TMLE_intercept <- private$TMLE.fit$TMLE_intercept

        if (!is.na(TMLE_intercept) && !is.nan(TMLE_intercept)) {
          update.Qstar.coef <- TMLE_intercept
        } else {
          update.Qstar.coef <- 0
        }

        # EIC_i_t_calc_unadjusted <- wts_TMLE * (Qkplus1 - Qk_hat)
        # print("EIC_i_t_calc_unadjusted"); print(mean(EIC_i_t_calc_unadjusted))

        ## Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        Qk_hat <- plogis(qlogis(Qk_hat) + update.Qstar.coef)
        iQ_all <- plogis(qlogis(iQ_all) + update.Qstar.coef)

        EIC_i_t_calc <- wts_TMLE * (Qkplus1 - Qk_hat)
        # print("EIC_i_t_calc_adjusted"); print(mean(EIC_i_t_calc))
        data$dat.sVar[self$idx_used_to_fit_initQ, ("EIC_i_t") := EIC_i_t_calc]
      }

      ## Q.k.hat is the prediction of the target parameter (\psi_hat) at the current time-point k (where we already set A(k) to A^*(k))
      ## This is either the initial G-COMP Q or the TMLE targeted version of the initial Q
      ## This prediction includes all newly censored observations.
      ## When stratifying Q fits, these predictions will also include all observations who just stopped following the treatment rule at time point k.
      data$dat.sVar[self$subset_idx, ("Qk_hat") := iQ_all]

      # Set the outcome for the next Q-regression: put Q[t] in (t-1).
      # only set the Qkplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
        ## using data.table:
        newsubset_idx <- (self$subset_idx - 1)

        ## The outcomes of these IDs (row numbers) need to be replaced with fold-specific outcomes for each training set
        ## Next training fold model can check if some of its indices overlap with 'train_idx_to_replace'
        ## How do we pass these predictions and store them?
        ## train_idx_to_replace <- (self$idx_used_to_fit_initQ - 1)

        data$dat.sVar[newsubset_idx, ("Qkplus1") := iQ_all]

        private$probA1 <- NULL
      } else {
        ## save prediction P(Q.kplus=1):
        private$probAeqa <- iQ_all
        private$probA1 <- NULL
      }

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      if (!self$keep_idx) self$wipe.all.indices # If we are planning on running iterative TMLE we will need the indicies used for fitting and predicting this Q
      # **********************************************************************
      invisible(self)
    },

    # **********************************************************************
    # Take a new TMLE fit and propagate it by first updating the Q(t) model and then updating the Q-model-based predictions for previous time-point
    # **********************************************************************
    Propagate_TMLE_fit = function(overwrite = TRUE, data, new.TMLE.fit, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitly asked
      TMLE_intercept <- new.TMLE.fit$TMLE_intercept
      if (!is.na(TMLE_intercept) && !is.nan(TMLE_intercept)) {
        update.Qstar.coef <- TMLE_intercept
      } else {
        update.Qstar.coef <- 0
      }

      ## Update the model predictions (Qk_hat) for initial Q[k] from GCOMP at time-point k.
      ## Based on TMLE update, the predictions now include ALL obs that are newly censored and just stopped following the rule at k:
      Qk_hat <- data$dat.sVar[self$subset_idx, "Qk_hat", with = FALSE][[1]]
      Qk_hat_updated <- plogis(qlogis(Qk_hat) + update.Qstar.coef)
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
      ## For CV-TMLE need to use holdout predictions
      ## These are also referred to as predictions from all validation splits, or holdouts or out-of-sample predictions
      # holdout <- self$CVTMLE
      assert_that(self$is.fitted)
      if (missing(newdata) && !is.null(private$probA1)) {
        ## probA1 will be a one column data.table, hence we extract and return the actual vector of predictions:
        return(private$probA1)
      } else if (missing(newdata) && is.null(private$probA1)) {
        probA1 <- gridisl::predict_SL(modelfit = private$model.fit,
                                      add_subject_data = FALSE,
                                      subset_idx = subset_idx,
                                      # use_best_retrained_model = TRUE,
                                      holdout = self$CVTMLE,
                                      verbose = gvars$verbose)

        ## probA1 will be a one column data.table, hence we extract and return the actual vector of predictions:
        private$probA1 <- probA1[[1]]
        return(private$probA1)

      } else {

        self$n <- newdata$nobs
        if (missing(subset_idx)) {
          subset_idx <- self$define.subset.idx(newdata, subset_exprs = self$subset_exprs)
        }

        if (!self$CVTMLE) {
          probA1 <- gridisl::predict_SL(modelfit = private$model.fit,
                                        newdata = newdata,
                                        add_subject_data = FALSE,
                                        subset_idx = subset_idx,
                                        # use_best_retrained_model = TRUE,
                                        holdout = FALSE,
                                        verbose = gvars$verbose)
        } else {
          probA1 <- gridisl::predict_SL(modelfit = private$model.fit,
                                        newdata = newdata,
                                        add_subject_data = FALSE,
                                        subset_idx = self$idx_used_to_fit_initQ,
                                        # use_best_retrained_model = TRUE,
                                        holdout = TRUE,
                                        verbose = gvars$verbose)
          probA1[, ("idx") := self$idx_used_to_fit_initQ]

          newObs_idx <- setdiff(subset_idx, self$idx_used_to_fit_initQ)
          if (length(newObs_idx)) {
            probA1_newObs <- gridisl::predict_SL(modelfit = private$model.fit,
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

        ## probA1 will be a one column data.table, hence we extract and return the actual vector of predictions:
        private$probA1 <- probA1[[1]]
        ## check that predictions P(A=1 | dmat) exist for all obs
        if (any(is.na(private$probA1) & !is.nan(private$probA1))) {
        # if (any(is.na(private$probA1) )) { # & !is.nan(private$probA1)
          stop("some of the predicted probabilities during seq Gcomp resulted in NAs, which indicates an error of a prediction routine")
        }

        ## Remove all modeling grid obj for xgboost (all the CV / training set models)
        ## This will still keep the best re-trained model object if method=="cv"/"holdout"
        ## Don't need to store these models, since we already made the prediction
        ## However, we will need these models if we want to report the model fit stats
        if (!(self$model_contrl[["fit_method"]] %in% "none")) private$model.fit$wipe.allmodels

        return(private$probA1)

      }
    },

    predictStatic = function(data, g0, gstar, subset_idx) {
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
      if (inherits(gcomp.pred.res, "try-error"))
        stop("error during prediction of the iterative GCOMP/TMLE")

      invisible(return(private$probA1))
    },

    predictStochastic = function(data, g0, gstar, subset_idx, stoch_indicator) {
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
        probA1 <- self$predictStatic(data, g0 = g0, gstar = gstar, subset_idx = subset_idx)

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
        probA1[subset_idx] <- probA1[subset_idx] * jointProb
        # Sum and keep looping
        stoch.probA1 <- stoch.probA1 + probA1[subset_idx]
      }
      stoch.probA1 <- stoch.probA1
      private$probA1[subset_idx] <- stoch.probA1
      return(invisible(private$probA1))
    },

    # Return the presaved prediction P(Q.kplus=1) only for a subset based on self$getsubset and private$probA1 (which has all n predictions)
    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      } else {
        stop("$predictAeqa should never be called for making predictions for GCOMP")
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
      # private$model.fit <- NULL
      # private$probAeqa <- NULL
      # self$binomialModelObj$emptydata
      # self$binomialModelObj$emptyY
      return(self)
    },
    wipe.all.indices = function() {
      self$idx_used_to_fit_initQ <- NULL
      self$subset_idx <- NULL
    },
    getfit = function() { private$model.fit },
    getTMLEfit = function() { private$TMLE.fit }
  ),

  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    TMLE.fit = list(NA),
    .outvar = NULL,
    probA1 = NULL,    # Predicted probA^s=1 conditional on Xmat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)