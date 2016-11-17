RegressionClassQlearn <- R6Class("RegressionClassQlearn",
  inherit = RegressionClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    Qreg_counter = integer(),
    t_period = integer(),
    TMLE = FALSE,
    stratifyQ_by_rule = FALSE,
    lower_bound_zero_Q = TRUE,
    skip_update_zero_Q = TRUE,
    regimen_names = NA,
    pool_regimes = FALSE,
    initialize = function(Qreg_counter,
                          t_period,
                          TMLE,
                          stratifyQ_by_rule,
                          regimen_names,
                          pool_regimes,
                          lower_bound_zero_Q = getopt("lower_bound_zero_Q"),
                          skip_update_zero_Q = getopt("skip_update_zero_Q"),
                          ...) {
      self$Qreg_counter <- Qreg_counter
      self$t_period <- t_period

      if (!missing(TMLE)) self$TMLE <- TMLE
      if (!missing(stratifyQ_by_rule)) self$stratifyQ_by_rule <- stratifyQ_by_rule
      if (!missing(regimen_names)) self$regimen_names <- regimen_names
      if (!missing(pool_regimes)) self$pool_regimes <- pool_regimes

      self$lower_bound_zero_Q <- lower_bound_zero_Q
      self$skip_update_zero_Q <- skip_update_zero_Q

      super$initialize(...)
    }
  ),
  active = list(
    get.reg = function() {
      list(Qreg_counter = self$Qreg_counter,
           t_period = self$t_period,
           TMLE = self$TMLE,
           outvar = self$outvar,
           predvars = self$predvars,
           outvar.class = self$outvar.class,
           subset_vars = self$subset_vars,
           subset_exprs = self$subset_exprs,
           subset_censored = self$subset_censored,
           stratifyQ_by_rule = self$stratifyQ_by_rule,
           lower_bound_zero_Q = self$lower_bound_zero_Q,
           regimen_names = self$regimen_names,
           pool_regimes = self$pool_regimes,
           model_contrl = self$model_contrl,
           censoring = self$censoring
           )
    }
  )
)

#************************************************
# TMLEs
#************************************************
tmle.update <- function(prev_Q.kplus1, init_Q_fitted_only, IPWts, lower_bound_zero_Q = TRUE, skip_update_zero_Q = TRUE) {
  QY.star <- NA
  if (sum(abs(IPWts)) < 10^-9) {
    update.Qstar.coef <- 0
    if (gvars$verbose) message("TMLE update cannot be performed since all IP-weights are exactly zero!")
    warning("TMLE update cannot be performed since all IP-weights are exactly zero!")
  } else if ((sum(prev_Q.kplus1[IPWts > 0]) < 10^-5) && skip_update_zero_Q) {
    update.Qstar.coef <- 0
  } else {
    #************************************************
    # TMLE update via weighted univariate ML (espsilon is intercept)
    #************************************************
    if (lower_bound_zero_Q) {
      prev_Q.kplus1[prev_Q.kplus1 < 10^(-4)] <- 10^(-4)
      init_Q_fitted_only[init_Q_fitted_only < 10^(-4)] <- 10^(-4)
    }
    off_TMLE <- qlogis(init_Q_fitted_only)

    m.Qstar <- try(speedglm::speedglm.wfit(X = matrix(1L, ncol=1, nrow=length(prev_Q.kplus1)),
                                          y = prev_Q.kplus1, weights = IPWts, offset = off_TMLE,
                                          # method=c('eigen','Cholesky','qr'),
                                          family = quasibinomial(), trace = FALSE, maxit = 1000),
                  silent = TRUE)

    if (inherits(m.Qstar, "try-error")) { # TMLE update failed
      if (gvars$verbose) message("attempt at running TMLE update with speedglm::speedglm.wfit has failed")
      warning("attempt at running TMLE update with speedglm::speedglm.wfit has failed")
      update.Qstar.coef <- 0
    } else {
      update.Qstar.coef <- m.Qstar$coef
    }
  }

  fit <- list(TMLE.intercept = update.Qstar.coef)
  class(fit)[2] <- "tmlefit"
  if (gvars$verbose) print("tmle update: " %+% update.Qstar.coef)
  return(fit)
}

## ---------------------------------------------------------------------
#' R6 Class for Q-Learning
#'
#'  R6 class for controlling the internal implementation of Q-learning functionality.
#'  Supports sequential (recursive) G-computation and longitudinal TMLE.
#'  Inherits from \code{BinaryOutcomeModel} R6 Class.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#'    \item{\code{regimen_names}} - \code{character}. Note used. For future pooling across regimens.
#'    \item{\code{classify}} - ... \code{logical}
#'    \item{\code{TMLE}} - ... \code{logical}
#'    \item{\code{nIDs}} - ... \code{integer}
#'    \item{\code{stratifyQ_by_rule}} - ... \code{logical}
#'    \item{\code{lower_bound_zero_Q}} - ... \code{logical}
#'    \item{\code{skip_update_zero_Q}} - ... \code{logical}
#'    \item{\code{Qreg_counter}} - ... \code{integer}
#'    \item{\code{t_period}} - ... \code{integer}
#'    \item{\code{idx_used_to_fit_initQ}} - ... \code{integer}
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, ...)}}{...}
#'   \item{\code{define.subset.idx(data, subset_vars, subset_exprs)}}{...}
#'   \item{\code{define_idx_to_fit_initQ(data)}}{...}
#'   \item{\code{define_idx_to_predictQ(data)}}{...}
#'   \item{\code{fit(overwrite = FALSE, data, ...)}}{...}
#'   \item{\code{Propagate_TMLE_fit(overwrite = TRUE, data, new.TMLE.fit, ...)}}{...}
#'   \item{\code{predict(newdata, subset_idx, ...)}}{...}
#'   \item{\code{predictStatic(data, g0, gstar, subset_idx)}}{...}
#'   \item{\code{predictStochastic(data, g0, gstar, subset_idx, stoch_indicator)}}{...}
#'   \item{\code{predictAeqa(newdata, ...)}}{...}
#'   \item{\code{get.fits()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{wipe.alldat}}{...}
#'    \item{\code{getfit}}{...}
#'    \item{\code{getTMLEfit}}{...}
#' }
#' @export
QlearnModel  <- R6Class(classname = "QlearnModel",
  inherit = BinaryOutcomeModel,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    regimen_names = character(), # for future pooling across regimens
    classify = FALSE,
    TMLE = TRUE,
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,
    lower_bound_zero_Q = TRUE,
    skip_update_zero_Q = TRUE,
    Qreg_counter = integer(), # Counter for the current sequential Q-regression (min is at 1)
    t_period = integer(),
    idx_used_to_fit_initQ = NULL,

    initialize = function(reg, ...) {
      super$initialize(reg, ...)
      self$Qreg_counter <- reg$Qreg_counter
      self$stratifyQ_by_rule <- reg$stratifyQ_by_rule
      self$lower_bound_zero_Q <- reg$lower_bound_zero_Q
      self$skip_update_zero_Q <- reg$skip_update_zero_Q
      self$t_period <- reg$t_period
      self$regimen_names <- reg$regimen_names
      self$TMLE <- reg$TMLE

      if (gvars$verbose) {print("initialized Q class"); reg$show()}

      invisible(self)
    },

    define.subset.idx = function(data, subset_vars, subset_exprs) {
      subset_idx <- data$evalsubst(subset_vars, subset_exprs)
      assert_that(is.logical(subset_idx))
      if ((length(subset_idx) < self$n) && (length(subset_idx) > 1L)) {
        if (gvars$verbose) message("subset_idx has smaller length than self$n; repeating subset_idx p times, for p: " %+% data$p)
        subset_idx <- rep.int(subset_idx, data$p)
        if (length(subset_idx) != self$n) stop("binomialModelObj$define.subset.idx: self$n is not equal to nobs*p!")
      }
      assert_that((length(subset_idx) == self$n) || (length(subset_idx) == 1L))
      return(subset_idx)
    },

    define_idx_to_fit_initQ = function(data) {
      # self$subset_vars is the name of the outcome var, checks for (!is.na(self$subset_vars))
      # self$subset_exprs is the time variable for selecting certain rows
      subset_idx <- self$define.subset.idx(data, subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)
      # excluded all censored observations:
      subset_idx <- subset_idx & data$uncensored_idx
      # if stratifying by rule, exclude all obs who are not following the rule:
      if (self$stratifyQ_by_rule) subset_idx <- subset_idx & data$rule_followers_idx
      return(subset_idx)
    },

    # Add all obs that were also censored at t and (possibly) those who just stopped following the rule
    define_idx_to_predictQ = function(data) {
      subset_idx <- self$define.subset.idx(data, subset_exprs = self$subset_exprs)
      return(subset_idx)
    },

    fit = function(overwrite = FALSE, data, iterTMLE = FALSE, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      # **********************************************************************
      # FITTING STEP OF Q-LEARNING
      # Select all obs at t who were uncensored & possibly were following the rule & had no missing outcomes (!is.na(self$subset_vars))
      # **********************************************************************
      self$subset_idx <- which(self$define_idx_to_fit_initQ(data))
      # save the subset used for fitting of the current initial Q[t] -> will be used for targeting
      self$idx_used_to_fit_initQ <- self$subset_idx
      # Fit model using Q.kplus as the outcome to obtain the inital model fit for Q[t]:
      private$model.fit <- self$binomialModelObj$fit(data, self$outvar, self$predvars, self$subset_idx, ...)
      self$is.fitted <- TRUE

      # **********************************************************************
      # PREDICTION STEP OF Q-LEARNING
      # Q prediction for everyone (including those who just got censored and those who just stopped following the rule)
      # **********************************************************************
      interventionNodes.g0 <- data$interventionNodes.g0
      interventionNodes.gstar <- data$interventionNodes.gstar
      # Determine which nodes are actually stochastic and need to be summed out
      stoch_indicator <- data$define.stoch.nodes(interventionNodes.gstar)
      any_stoch <- sum(stoch_indicator) > 0

      # For prediction need to add all obs that were also censored at t and (possibly) those who just stopped following the rule
      self$subset_idx <- which(self$define_idx_to_predictQ(data))

      if (!any_stoch) {
        probA1 <- self$predictStatic(data, g0 = interventionNodes.g0,
                                           gstar = interventionNodes.gstar,
                                           subset_idx = self$subset_idx)
      } else {
        # For all stochastic nodes, need to integrate out w.r.t. the support of each node
        probA1 <- self$predictStochastic(data, g0 = interventionNodes.g0,
                                               gstar = interventionNodes.gstar,
                                               subset_idx = self$subset_idx,
                                               stoch_indicator = stoch_indicator)
      }
      init_Q_all_obs <- probA1[self$subset_idx]

      # ------------------------------------------------------------------------------------------------------------------------
      # Alternative to above:
      # ------------------------------------------------------------------------------------------------------------------------
      # 1. reset current exposures at t to their counterfactual values (e.g., by renaming the columns)
      # note: this step is unnecessary when doing fitting only among rule-followers
      # self$binomialModelObj$replaceCols(data, self$subset_idx, data$nodes$, , ...)
      # 2. predict for all obs who were used for fitting t (re-using the same design mat used for fitting)
      # probAeqa <- self$predictAeqa(...)
      # self$predict(...)
      # 3. add observations that were censored at current t to the subset and do predictions for them
      # ------------------------------------------------------------------------------------------------------------------------
      # *** NEED TO ADD: possibly intergrate out over the support of the stochastic intervention (weighted sum of P(A(t)=a(t)))
      # ------------------------------------------------------------------------------------------------------------------------
      # *** NEED TO ADD: do MC sampling to perform the same integration for stochastic g^*
      # ------------------------------------------------------------------------------------------------------------------------
      # print("initial mean(Q.kplus1) for ALL obs at t=" %+% self$t_period %+% ": " %+% round(mean(init_Q_all_obs), 4))

      # **********************************************************************
      # TARGETING STEP OF THE TMLE
      # the TMLE update is performed only among obs who were involved in fitting of the initial Q above (self$idx_used_to_fit_initQ)
      # **********************************************************************
      # Predicted outcome from the previous Seq-GCOMP/TMLE iteration, was saved in the current row t
      # If the person just failed at t this is always 1 (deterministic)
      prev_Q.kplus1 <- data$dat.sVar[self$idx_used_to_fit_initQ, "Q.kplus1", with = FALSE][[1]]
      data$dat.sVar[self$idx_used_to_fit_initQ, "prev_Q.kplus1" := eval(as.name("Q.kplus1"))]
      # data$dat.sVar[self$idx_used_to_fit_initQ, ]

      if (self$TMLE) {
        # TMLE offset will be based on the initial prediction of Q among fitted only (log(x/[1-x])):
        init_Q_fitted_only <- probA1[self$idx_used_to_fit_initQ]

        # tmle update will error out if some predictions are exactly 0:
        off_TMLE <- qlogis(init_Q_fitted_only)
        # print("initial mean(Q.kplus1) among fitted obs only at t=" %+% self$t_period %+% ": " %+% round(mean(init_Q_all_obs), 4))
        # Cumulative IPWeights for current t:
        wts_TMLE <- data$IPwts_by_regimen[self$idx_used_to_fit_initQ, "cum.IPAW", with = FALSE][[1]]

        # TMLE update based on the IPWeighted logistic regression model with offset and intercept only:
        private$TMLE.fit <- tmle.update(prev_Q.kplus1 = prev_Q.kplus1,
                                        init_Q_fitted_only = init_Q_fitted_only,
                                        IPWts = wts_TMLE,
                                        lower_bound_zero_Q = self$lower_bound_zero_Q,
                                        skip_update_zero_Q = self$skip_update_zero_Q)
        TMLE.intercept <- private$TMLE.fit$TMLE.intercept
        # TMLE.cleverCov.coef <- private$TMLE.fit$TMLE.cleverCov.coef
        # print("TMLE Intercept: " %+% round(TMLE.intercept, 5))
        # print("TMLE Intercept: " %+% TMLE.intercept)

        if (!is.na(TMLE.intercept) && !is.nan(TMLE.intercept)) {
          update.Qstar.coef <- TMLE.intercept
        } else {
          update.Qstar.coef <- 0
        }

        # Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        init_Q_fitted_only <- plogis(qlogis(init_Q_fitted_only) + update.Qstar.coef)
        init_Q_all_obs <- plogis(qlogis(init_Q_all_obs) + update.Qstar.coef)
        # wts_TMLE_all_obs <- data$IPwts_by_regimen[self$subset_idx, "cum.IPAW", with = FALSE][[1]]
        # init_Q_all_obs <- plogis(qlogis(init_Q_all_obs) + TMLE.cleverCov.coef * wts_TMLE_all_obs)
        # print("TMLE update of mean(Q.kplus1) at t=" %+% self$t_period %+% ": " %+% mean(init_Q_all_obs))]

        EIC_i_t_calc <- wts_TMLE * (prev_Q.kplus1 - init_Q_fitted_only)
        data$dat.sVar[self$idx_used_to_fit_initQ, ("EIC_i_t") := EIC_i_t_calc]
      }

      # Save all predicted vals as Q.kplus1[t] in row t or first target and then save targeted values:
      # rowidx_t <- which(self$subset_idx)
      data$dat.sVar[self$subset_idx, "Q.kplus1" := init_Q_all_obs]

      # Set the outcome for the next Q-regression: put Q[t] in (t-1), this will be overwritten with next prediction
      # only set the Q.kplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
      # rowidx_t.minus1 <- self$subset_idx - 1
        data$dat.sVar[(self$subset_idx - 1), "Q.kplus1" := init_Q_all_obs]
        private$probA1 <- NULL
      } else {
        # save prediction P(Q.kplus=1) only for a subset in self$getsubset that was used for prediction
        private$probAeqa <- init_Q_all_obs
        private$probA1 <- NULL
      }

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      if (!iterTMLE) self$wipe.all.indices # If we are planning on running iterative TMLE we will need the indicies used for fitting and predicting this Q
      # **********************************************************************
      invisible(self)
    },

    # **********************************************************************
    # Take a new TMLE fit and propagate it by first updating the Q(t) model and then updating the Q-model-based predictions for previous time-point
    # **********************************************************************
    Propagate_TMLE_fit = function(overwrite = TRUE, data, new.TMLE.fit, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked
      TMLE.intercept <- new.TMLE.fit$TMLE.intercept
      if (!is.na(TMLE.intercept) && !is.nan(TMLE.intercept)) {
        update.Qstar.coef <- TMLE.intercept
      } else {
        update.Qstar.coef <- 0
      }
      # rowidx_t <- which(self$subset_idx)
      # Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
      Q.kplus1 <- data$dat.sVar[self$subset_idx, "Q.kplus1", with = FALSE][[1]]
      Q.kplus1.new <- plogis(qlogis(Q.kplus1) + update.Qstar.coef)
      # Save all predicted vals as Q.kplus1[t] in row t or first target and then save targeted values:
      data$dat.sVar[self$subset_idx, "Q.kplus1" := Q.kplus1.new]
      # Set the outcome for the next Q-regression: put Q[t] in (t-1), this will be overwritten with next prediction
      # only set the Q.kplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
        # rowidx_t.minus1 <- rowidx_t - 1
        data$dat.sVar[(self$subset_idx-1), "prev_Q.kplus1" := Q.kplus1.new]
        private$probA1 <- NULL
      } else {
        # save prediction P(Q.kplus=1) only for a subset in self$getsubset that was used for prediction
        private$probAeqa <- Q.kplus1.new
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
      assert_that(self$is.fitted)
      if (missing(newdata) && !is.null(private$probA1)) {
        return(private$probA1)
      } else if (missing(newdata) && is.null(private$probA1)) {
        private$probA1 <- self$binomialModelObj$predictP1(subset_idx = subset_idx)
        return(private$probA1)
      } else {
        self$n <- newdata$nobs
        if (missing(subset_idx)) {
          subset_idx <- self$define.subset.idx(newdata, subset_exprs = self$subset_exprs)
        }
        private$probA1 <- self$binomialModelObj$predictP1(data = newdata, subset_idx = subset_idx)
        if (any(is.na(private$probA1[subset_idx]) & !is.nan(private$probA1[subset_idx]))) {
          stop("some of the predicted probabilities during seq Gcomp resulted in NAs, which indicates an error of a prediction routine")
        }
        # assert_that(!any(is.na(private$probA1[subset_idx]))) # check that predictions P(A=1 | dmat) exist for all obs.

        invisible(return(private$probA1))
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

      if (inherits(gcomp.pred.res, "try-error")) { # prediction in seq-Gcomp has failed
        stop("attempt at prediction during GCOMP/TMLE has failed")
      }
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
      invisible(return(private$probA1))
    },

    # Return the presaved prediction P(Q.kplus=1) only for a subset based on self$getsubset and private$probA1 (which has all n predictions)
    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      } else {
        stop("$predictAeqa should never be called for making actual predictions")
      }
      # if (is.null(self$getsubset)) stop("cannot make predictions after self$subset_idx is erased (set to NULL)")
      # invisible(return(private$probAeqa))
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
      # private$probA1 <- NULL
      # private$probAeqa <- NULL

      self$binomialModelObj$emptydata
      self$binomialModelObj$emptyY
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