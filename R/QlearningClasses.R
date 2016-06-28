RegressionClassQlearn <- R6Class("RegressionClassQlearn",
  inherit = RegressionClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    Qreg_counter = integer(),
    t_period = integer(),
    TMLE = FALSE,
    regimen_names = NA,
    stratifyQ_by_rule = FALSE,
    pool_regimes = FALSE,
    initialize = function(Qreg_counter, t_period, TMLE, regimen_names, stratifyQ_by_rule, pool_regimes, ...) {
      self$Qreg_counter <- Qreg_counter
      self$t_period <- t_period

      if (!missing(TMLE)) self$TMLE <- TMLE
      if (!missing(regimen_names)) self$regimen_names <- regimen_names
      if (!missing(stratifyQ_by_rule)) self$stratifyQ_by_rule <- stratifyQ_by_rule
      if (!missing(pool_regimes)) self$pool_regimes <- pool_regimes
      super$initialize(...)
    }
  ),
  active = list(
    get.reg = function() {
      list(Qreg_counter = self$Qreg_counter,
           t_period = self$t_period,
           TMLE = self$TMLE,
           regimen_names = self$regimen_names,
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

#************************************************
# TMLEs
#************************************************
tmle.update <- function(prev_Q.kplus1, off, IPWts, determ.Q, predictQ = TRUE) {
  QY.star <- NA
  #************************************************
  # TMLE update via weighted univariate ML (espsilon is intercept)
  #************************************************
  # SuppressGivenWarnings(
  update_t <- system.time(
    m.Qstar <- speedglm::speedglm.wfit(X = matrix(1L, ncol=1, nrow=length(prev_Q.kplus1)),
                                        y = prev_Q.kplus1, weights = IPWts, offset = off,
                                        family = quasibinomial(), trace = FALSE, maxit = 1000)    # ,
    )
  print("time to perform tmle update:"); print(update_t)
    # GetWarningsToSuppress(TRUE))
  update.Qstar.coef <- m.Qstar$coef
  # m.Qstar.fit <- list(coef = m.Qstar$coef, linkfun = "logit_linkinv", fitfunname = "speedglm")
  # class(m.Qstar.fit) <- c(class(m.Qstar.fit), c("speedglmS3"))
  # SuppressGivenWarnings(m.Qstar <- glm(Y ~ offset(off), data = data.frame(Y = Y, off = off), weights = IPWts,
                                            # subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))

  # update.Qstar.coef <- coef(m.Qstar)
  # if (predictQ) {
  #   QY.star <- prev_Q.kplus1
  #   if (!is.na(update.Qstar.coef)) QY.star <- plogis(off + update.Qstar.coef)
  # }

  #************************************************
  # (DISABLED) g_IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
  # 02/16/13: IPTW estimator (Y_i * prod_{j \\in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
  # g_wts <- iptw_est(k = est_params_list$Kmax, data = data, node_l = nodes, m.gN = est_params_list$m.g0N,
  #                      f.gstar = est_params_list$f.gstar, f.g_args = est_params_list$f.g_args, family = "binomial",
  #                      NetInd_k = est_params_list$NetInd_k, lbound = est_params_list$lbound, max_npwt = est_params_list$max_npwt,
  #                      f.g0 = est_params_list$f.g0, f.g0_args = est_params_list$f.g0_args)
  # Y_IPTW_g <- Y
  # Y_IPTW_g[!determ.Q] <- Y[!determ.Q] * g_wts[!determ.Q]
  return(list(update.Qstar.coef = update.Qstar.coef, QY.star = QY.star))
}

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
    Qreg_counter = integer(), # Counter for the current sequential Q-regression (min is at 1)
    t_period = integer(),
    idx_used_to_fit_initQ = NULL,

    initialize = function(reg, ...) {
      super$initialize(reg, ...)
      self$stratifyQ_by_rule <- reg$stratifyQ_by_rule
      self$Qreg_counter <- reg$Qreg_counter
      self$t_period <- reg$t_period
      self$regimen_names <- reg$regimen_names
      self$TMLE <- reg$TMLE
      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and save those to input data
    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$nIDs <- data$nuniqueIDs

      if (gvars$verbose) print("fitting G-COMP model: " %+% self$show())
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      # **********************************************************************
      # FITTING STEP OF Q-LEARNING
      # **********************************************************************
      # select all obs at current t, with no missing outcome (!is.na(self$subset_vars))
      self$subset_idx <- self$define.subset.idx(data, subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)
      # excluded all censored observations:
      self$subset_idx <- self$subset_idx & data$uncensored_idx
      # if stratifying by rule, exclude all obs who are not following the rule:
      if (self$stratifyQ_by_rule) self$subset_idx <- self$subset_idx & data$rule_followers_idx
      # save the subset used for fitting of the current initial Q[t] -> will be used for targeting
      self$idx_used_to_fit_initQ <- self$subset_idx
      # obtain the inital model fit for Q[t]:
      private$model.fit <- self$binomialModelObj$fit(data, self$outvar, self$predvars, self$subset_idx, ...)
      self$is.fitted <- TRUE
      print("Q-learning for: " %+% self$subset_exprs); print(self$get.fits())

      # **********************************************************************
      # PREDICTION STEP OF Q-LEARNING
      # Q prediction for everyone (including those who just got censored and those who just stopped following the rule)
      # **********************************************************************
      # Set current A's to the counterfactual exposures in the data (for predicting Q):
      interventionNodes.g0 <- data$interventionNodes.g0
      interventionNodes.gstar <- data$interventionNodes.gstar
      data$swapNodes(current = interventionNodes.g0, target = interventionNodes.gstar)

      # Add to design matrix all obs that were also censored at t and (possibly) those just stopped following the rule:
      self$subset_idx <- self$define.subset.idx(data, subset_exprs = self$subset_exprs)
      print("performing initial Q-prediction for N = " %+% sum(self$subset_idx))

      # Predict Prob(Q.init = 1) for all observations in subset_idx (note: probAeqa is never used, only private$probA1)
      probAeqa <- self$predictAeqa(data, subset_idx = self$subset_idx)

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

      # Obtain the initial Q prediction P(Q=1|...) for EVERYBODY (including those who just got censored and those who just stopped following the rule)
      init_Q_all_obs <- private$probA1[self$subset_idx]
      print("initial mean(Q.kplus1) for ALL obs at t=" %+% self$t_period %+% ": " %+% round(mean(init_Q_all_obs), 4))
      # **********************************************************************
      # TARGETING STEP OF THE TMLE
      # the TMLE update is performed only among obs who were involved in fitting of the initial Q above (self$idx_used_to_fit_initQ)
      # **********************************************************************
      if (self$TMLE) {
        # Predicted outcome from the previous Seq-GCOMP/TMLE iteration, was saved in the current row t
        # If the person just failed at t this is always 1 (deterministic)
        prev_Q.kplus1 <- data$dat.sVar[self$idx_used_to_fit_initQ, "Q.kplus1", with = FALSE][[1]]
        # TMLE offset based will be based on the initial prediction of Q above log(x/[1-x]):
        init_Q_fitted_only <- private$probA1[self$idx_used_to_fit_initQ]
        off_TMLE <- qlogis(init_Q_fitted_only)
        print("initial mean(Q.kplus1) among fitted obs only at t=" %+% self$t_period %+% ": " %+% round(mean(init_Q_all_obs), 4))

        # Cumulative IPWeights for current t:
        wts_TMLE <- data$IPwts_by_regimen[self$idx_used_to_fit_initQ, "cumm.IPAW", with = FALSE][[1]]

# summary(off_TMLE)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -36.0400 -22.1400  -4.8970 -12.4900  -3.5830   0.4458
        # summary(init_Q_fitted_only)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.000000 0.007413 0.016660 0.027030 0.609600

        # TMLE update based on the IPWeighted logistic regression model with offset and intercept only:
        tmle_res <- tmle.update(prev_Q.kplus1 = prev_Q.kplus1, off = off_TMLE, IPWts = wts_TMLE)
        print("TMLE Intercept: " %+% round(tmle_res$update.Qstar.coef, 5))
        if (!is.na(tmle_res$update.Qstar.coef)) {
          update.Qstar.coef <- tmle_res$update.Qstar.coef
        } else {
          update.Qstar.coef <- 0
        }

        # Updated the model predictions (Q.star) for init_Q based on TMLE update using ALL obs (inc. newly censored and newly non-followers):
        init_Q_all_obs <- plogis(qlogis(init_Q_all_obs) + tmle_res$update.Qstar.coef)
        # Q.kplus1 <- plogis(off_TMLE + update.Qstar.coef)
        print("TMLE update of mean(Q.kplus1) at t=" %+% self$t_period %+% ": " %+% mean(init_Q_all_obs))
      }

      # Save all predicted vals as Q.kplus1[t] in row t or first target and then save targeted values:
      rowidx_t <- which(self$subset_idx)
      data$dat.sVar[rowidx_t, "Q.kplus1" := init_Q_all_obs]

      # Set the outcome for the next Q-regression: put Q[t] in (t-1), this will be overwritten with next prediction
      # only set the Q.kplus1 while self$Qreg_counter > 1, self$Qreg_counter == 1 implies that Q-learning finished & reached the minimum/first time-point period
      if (self$Qreg_counter > 1) {
        rowidx_t.minus1 <- rowidx_t - 1
        data$dat.sVar[rowidx_t.minus1, "Q.kplus1" := init_Q_all_obs]
        private$probA1 <- NULL
        private$probAeqa <- NULL
      }

      # Reset back the observed exposure to A[t] (swap back by renaming columns)
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
      self$idx_used_to_fit_initQ <- NULL
      self$subset_idx <- NULL
      self$binomialModelObj$emptydata
      self$binomialModelObj$emptyY
      return(self)
    }
  )
)