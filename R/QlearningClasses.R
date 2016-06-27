RegressionClassQlearn <- R6Class("RegressionClassQlearn",
  inherit = RegressionClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    Qreg_counter = integer(),
    t_period = integer(),
    stratifyQ_by_rule = FALSE,
    pool_regimes = FALSE,
    initialize = function(Qreg_counter, t_period, stratifyQ_by_rule, pool_regimes, ...) {
      self$Qreg_counter <- Qreg_counter
      self$t_period <- t_period
      if (!missing(stratifyQ_by_rule)) self$stratifyQ_by_rule <- stratifyQ_by_rule
      if (!missing(pool_regimes)) self$pool_regimes <- pool_regimes
      super$initialize(...)
    }
  ),
  active = list(
    get.reg = function() {
      list(Qreg_counter = self$Qreg_counter,
           t_period = self$t_period,
           outvar = self$outvar,
           predvars = self$predvars,
           outvar.class = self$outvar.class,
           subset_vars = self$subset_vars,
           subset_exprs = self$subset_exprs,
           subset_censored = self$subset_censored,
           stratifyQ_by_rule = self$stratifyQ_by_rule,
           pool_regimes = self$pool_regimes,
           model_contrl = self$model_contrl,
           censoring = self$censoring
           )
    }
  )
)

#' @export
QlearnModel  <- R6Class(classname = "QlearnModel",
  inherit = BinaryOutcomeModel,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    classify = FALSE,
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,
    Qreg_counter = integer(), # Counter for the current sequential Q-regression (min is at 1)
    t_period = integer(),

    initialize = function(reg, ...) {
      super$initialize(reg, ...)
      self$stratifyQ_by_rule <- reg$stratifyQ_by_rule
      self$Qreg_counter <- reg$Qreg_counter
      self$t_period <- reg$t_period
      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and save those to input data
    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
    # , predict = FALSE,
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs

      if (gvars$verbose) print("fitting G-COMP the model: " %+% self$show())
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      # **********************************************************************
      # FITTING STEP OF Q-LEARNING
      # **********************************************************************
      # select all observations for current t==self$t_period and possibly (!self$gstar > 0)

      # select all obs at current t, with no missing outcome (!is.na(self$subset_vars))
      self$subset_idx <- self$define.subset.idx(data, subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)

      # excluded all censored observations:
      self$subset_idx <- self$subset_idx & data$uncensored_idx

      # if stratifying by rule, exclude all obs who are not following the rule:
      if (self$stratifyQ_by_rule) self$subset_idx <- self$subset_idx & data$rule_followers_idx

      private$model.fit <- self$binomialModelObj$fit(data, self$outvar, self$predvars, self$subset_idx, ...)

      self$is.fitted <- TRUE
      print("Q-learning for: " %+% self$subset_exprs); print(self$get.fits())

      # **********************************************************************
      # PREDICTION STEP OF Q-LEARNING
      # **********************************************************************
      # 1.
      # Set A's to counterfactuals in data, then reset the whole design matrix with a new subset that includes everyone who just censored at t.
      interventionNodes.g0 <- data$interventionNodes.g0
      interventionNodes.gstar <- data$interventionNodes.gstar
      data$swapNodes(current = interventionNodes.g0, target = interventionNodes.gstar)
      self$subset_idx <- self$define.subset.idx(data, subset_exprs = self$subset_exprs)

      print("Q-prediction for N: " %+% sum(self$subset_idx))
      probAeqa <- self$predictAeqa(data, subset_idx = self$subset_idx)

      # Alternative to 1:
      # *** NEED TO: 1. reset current exposures at t to their counterfactual values (e.g., by renaming the columns)
      # note: this step is unnecessary when doing fitting only among rule-followers
      # self$binomialModelObj$replaceCols(data, self$subset_idx, data$nodes$, , ...)
      # 2. predict for all obs who were used for fitting t (re-using the same design mat used for fitting)
      # probAeqa <- self$predictAeqa(...)
      # self$predict(...)
      # *** NEED TO: 3. add observations that were censored at current t to the subset and do predictions for them

      # ------------------------------------------------------------------------------------------------------------------------
      # *** NEED TO ADD: possibly intergrate out over the support of the stochastic intervention (weighted sum of P(A(t)=a(t)))
      # ------------------------------------------------------------------------------------------------------------------------
      # *** NEED TO ADD: do MC sampling to perform the same integration
      # ------------------------------------------------------------------------------------------------------------------------

      # 2. save all predicted vals as Q.kplus1[t] in row t (saved for later targeting/etc):
      rowidx_t <- which(self$subset_idx)
      data$dat.sVar[rowidx_t, "Q.kplus1" := private$probA1[self$subset_idx]]
      # data$dat.sVar[1:40,]

      # 3. set the outcome for the next Q-regression: put Q[t] in (t-1), this will be overwritten with next prediction
      #    only set the Q.kplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
        rowidx_t.minus1 <- rowidx_t - 1
        data$dat.sVar[rowidx_t.minus1, "Q.kplus1" := private$probA1[self$subset_idx]]
        private$probA1 <- NULL
        private$probAeqa <- NULL
      }

      # 4. reset back the observed exposure to A[t] (swap back by renaming columns)
      data$swapNodes(current = interventionNodes.gstar, target = interventionNodes.g0)

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$model.fit to generate predictions for data:
    predict = function(newdata, subset_idx, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata) && !is.null(private$probA1)) {
        return(private$probA1)
      } else if (missing(newdata) && is.null(private$probA1)) {
        private$probA1 <- self$binomialModelObj$predictP1(subset_idx = self$subset_idx)
        return(private$probA1)
      } else {
        self$n <- newdata$nobs
        if (missing(subset_idx)) {
          self$subset_idx <- self$define.subset.idx(newdata, subset_exprs = self$subset_exprs)
        } else {
          self$subset_idx <- subset_idx
        }
        private$probA1 <- self$binomialModelObj$predictP1(data = newdata, subset_idx = self$subset_idx)
        return(private$probA1)
      }
    },

    # Predict P(Q=1)
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, subset_idx, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      }
      if (is.null(self$getsubset)) stop("cannot make predictions after self$subset_idx is erased (set to NULL)")
      self$predict(newdata = newdata, subset_idx = subset_idx)
      # probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      probA1 <- private$probA1[self$getsubset]
      assert_that(!any(is.na(probA1))) # check that predictions P(A=1 | dmat) exist for all obs.
      # probAeqa[self$getsubset] <- probA1
      private$probAeqa <- probA1
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # self$wipe.alldat
      # **********************************************************************
      return(probA1)
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

    # Returns the object that contains the actual model fits (itself)
    get.fits = function() {
      model.fit <- self$getfit
      return(list(model.fit))
    }
  ),

  active = list(
    wipe.alldat = function() {
      # private$probA1 <- NULL
      # private$probAeqa <- NULL
      self$subset_idx <- NULL
      self$binomialModelObj$emptydata
      self$binomialModelObj$emptyY
      return(self)
    }
  )
)