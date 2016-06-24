# SGCompRegClass <- R6Class("SGCompRegClass",
#   inherit = RegressionClass,
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     # outvar = character(),          # the outcome variable name (Ynode)
#     # outvar.class = character(),
#     Qforms = character(),
#     # n_regs = integer(),
#     # n_timepts = integer(),       # number of time points
#     # nodes = list(),
#     # predvars = list(),             # list of predictors for each regression from Qforms
#     # not used for now:
#     Anodes = list(),               # list of treatment var names, by timepoint
#     Cnodes = list(),               # list of censoring var names, by timepoint
#     Lnodes = list(),               # list of time-varying confounder names, by timepoint
#     Wnodes = character(),          # character vector of baseline covariate names
#     subset = NULL,                 # subset expression (later evaluated to logical vector in the envir of the data)
#     ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
#     pool = logical()
#   )
# )

RegressionClassQlearn <- R6Class("RegressionClassQlearn",
  inherit = RegressionClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    Qreg_counter = integer(),
    t_period = integer(),
    # subset_censored = logical(),
    stratifyQ_by_rule = FALSE,
    pool_regimes = FALSE,
    initialize = function(Qreg_counter, t_period, stratifyQ_by_rule, pool_regimes, ...) {
      self$Qreg_counter <- Qreg_counter
      self$t_period <- t_period
      # subset_censored
      # if (!missing(subset_censored)) self$subset_censored <- subset_censored
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
    nIDs = integer(),
    stratifyQ_by_rule = FALSE,
    # interventionNodes.g0 = character(),
    # interventionNodes.gstar = character(),
    # subset_t_expr = NULL,
    # subset_rule_expr = NULL,
    # subset_CENS_expr = NULL,  # THE LOGICAL EXPRESSION, WHEN EVALUTED IDENTIFIES ALL CENSORED OBS (at current t)
    # subset_t_idx = NULL,
    # subset_rule_idx = NULL,
    # subset_CENS_idx = NULL,

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
      private$model.fit$params <- self$show(print_format = FALSE)
      self$is.fitted <- TRUE

      print("Q-learning for: " %+% self$subset_exprs); print(self$get.fits())

      # browser()

      # **********************************************************************
      # PREDICTION STEP OF Q-LEARNING
      # **********************************************************************
      # 1.
      # Set A's to counterfactuals in data, then reset the whole design matrix with a new subset that includes everyone who just censored at t.
      interventionNodes.g0 <- data$interventionNodes.g0
      interventionNodes.gstar <- data$interventionNodes.gstar
      data$swapNodes(current = interventionNodes.g0, target = interventionNodes.gstar)
      data$dat.sVar
      self$subset_idx <- self$define.subset.idx(data, subset_exprs = self$subset_exprs)
      probAeqa <- self$predictAeqa(data, subset_idx = self$subset_idx)

      # Alternative 1b:
      # *** NEED TO: 1. reset current exposures at t to their counterfactual values (e.g., by renaming the columns)
      # note: this step is unnecessary when doing fitting only among rule-followers
      # self$binomialModelObj$replaceCols(data, self$subset_idx, data$nodes$, , ...)
      # 2. predict for all obs who were used for fitting t (re-using the same design mat used for fitting)
      # probAeqa <- self$predictAeqa(...)
      # self$predict(...)
      # *** NEED TO: 3. add observations that were censored at current t to the subset and do predictions for them




      # *** NEED TO: 2a. possibly intergrate out over the support of the stochastic intervention (weighted sum of P(A(t)=a(t)))
      # *** NEED TO: 2b. do MC sampling to perform the same integration
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

      # *** NEED TO: 4. reset back the observed exposure to A[t] (renaming columns)
      # ...

      # sum(self$subset_idx)
      # self$t_period
      # self$Qreg_counter
      # data$dat.sVar[1:50, ]
      # browser()

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

#' @export
fitSeqGcomp <- function(OData,
                        t = OData$max.t,
                        Qforms,
                        gstar_TRT = NULL,
                        gstar_MONITOR = NULL,
                        stratifyQ_by_rule = FALSE,
                        rule_followers_colname = NULL,
                        params = list(),
                        verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  assert_that(is.list(params))


  # ------------------------------------------------------------------------------------------------
  # TO DO:
  # ------------------------------------------------------------------------------------------------
  # *) Get exactly the same est stratifyQ_by_rule TRUE or FALSE -> MUST BE AN ERROR
  # *) Need to redefine gstar_TRT into a vector of counterfactual probabilities for TRT nodes
  # *) gstar_MONITOR is already defined correctly, need to allow it to be multivariate (more than one node)
  # *) allow each gstar to be vector (like abar=(0,0,0,0)) or a matrix of counterfactual treatments
  # *) finally, when <1 and >0, assume its the counterfactual probability P(TRT^*=1)=gstar_TRT
  # *) rule_followers_colname: RULE FOLLOWERS COLUMN NEEDS TO BE EVALUTED AUTOMATICALLY!!!!
  # *) The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------------------------
  # **** Define the intervention nodes
  # ------------------------------------------------------------------------------------------------
  if (!is.null(gstar_TRT)) {
    gstar.A <- gstar_TRT
  } else {
    gstar.A <- nodes$Anodes # use the actual observed exposure (no intervention on TRT)
  }
  if (!is.null(gstar_MONITOR)) {
    gstar.N <- gstar_MONITOR
  } else {
    gstar.N <- nodes$Nnodes # use the actual observed monitoring probability (no intervention on MONITOR)
  }

  interventionNodes.g0 <- c(nodes$Anodes, nodes$Nnodes)
  interventionNodes.gstar <- c(gstar.A, gstar.N)
  OData$interventionNodes.g0 <- interventionNodes.g0
  OData$interventionNodes.gstar <- interventionNodes.gstar

  # ------------------------------------------------------------------------------------------------
  # **** Set the uncensored and rule nodes
  # ------------------------------------------------------------------------------------------------
  OData$uncensored_idx <- OData$dat.sVar[, list(uncensored_idx = as.logical(rowSums(.SD, na.rm = TRUE) == eval(OData$noCENScat))), .SDcols = nodes$Cnodes][["uncensored_idx"]]
  # !!!! NOTE THIS LINE IS TECHNICALLY INCORRECT: RULE FOLLOWERS COLUMN NEEDS TO BE EVALUTED AUTOMATICALLY!!!!
  if (stratifyQ_by_rule) {
    assert_that(!is.null(rule_followers_colname))
    OData$rule_followers_idx <- OData$dat.sVar[, list(rule_followers_idx = .SD > 0), .SDcols = rule_followers_colname][["rule_followers_idx"]]
  }

  # ------------------------------------------------------------------------------------------------
  # **** TO DO: The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------
  Qperiods <- rev(OData$min.t:t)
  Qreg_idx <- rev(seq_along(Qperiods))
  stratify_Q <- as.list(nodes[['tnode']] %+% " == " %+% (Qperiods))
  names(stratify_Q) <- rep.int("Q.kplus1", length(stratify_Q))

  # ------------------------------------------------------------------------------------------------
  # Process the input formulas and stratification settings;
  # Define regression classes for Q.Y and put them in a single list of regressions.
  # ------------------------------------------------------------------------------------------------
  Qforms.default <- rep.int("Q.kplus1 ~ Lnodes + Anodes + Cnodes + Nnodes", length(Qperiods))

  # ------------------------------------------------------------------------------------------------
  # TMLE:
  # Initiate the Q.kplus1 - need to do this for each regimen
  # That column keeps the tabs on the running Q-fit (SEQ G-COMP)
  OData$dat.sVar[, "Q.kplus1" := as.numeric(get(OData$nodes$Ynode))]
  OData$def.types.sVar()
  # ------------------------------------------------------------------------------------------------
  if (missing(Qforms)) Qforms <- Qforms.default

  Q_regs_list <- vector(mode = "list", length = length(stratify_Q))
  names(Q_regs_list) <- unlist(stratify_Q)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")

  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    reg <- RegressionClassQlearn$new(Qreg_counter = Qreg_idx[i], t_period = Qperiods[i],
                                     outvar = "Q.kplus1", predvars = regform$predvars, outvar.class = list("Qlearn"),
                                     subset_vars = list("Q.kplus1"), subset_exprs = stratify_Q[i], model_contrl = params,
                                     censoring = FALSE)
    Q_regs_list[[i]] <- reg
  }

  # browser()

  Qlearn.fit <- GenericModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)

  # ------------------------------------------------------------------------------------------
  # Outstanding issues:
  # ------------------------------------------------------------------------------------------
  # *** Need to define A^* at the begining, as a column (one for each regimen) -> should be returned by defineTRTrules()
  # *** Censoring at t: Should be flexible and allowed to be defined as either those who are censored at t or those who are no longer following the rule at t.
  # *** Need to add to current subsets expression for conditioning on uncensored observation (or censored) at each t or rule non-followers
  #     Q_regs_list <- stratify_by_uncensored(Q_regs_list)
  #     OData$dat.sVar[t == 0L, ]
  # ------------------------------------------------------------------------------------------
  # *** Stratification/pooling ***
  #     Since new regimen results in new Q.kplus1 and hence new outcome -> requires a separate regression for each regimen
  #     => can either pool all Q.regimens at once (fit one model for repated observations, one for each rule, smoothing over rules)
  #     => can use the same stacked dataset of regime-followers, but with a separate stratification for each regimen.
  # ------------------------------------------------------------------------------------------
  # *** Accessing correct QlearnModel ***
  #     Need to come up with a good way of accessing correct terminal QlearnModel to get final Q-predictions in private$probAeqa
  #     These predictions are for the final n observations that were used for evaluating E[Y^a]
  #     However, if we do stratify=TRUE and use a stack version of g-comp with different strata all in one stacked database it will be difficult
  #     to pull the right QlearnModel.
  #     Alaternative is to ignore that completely and go directly for the observed data (by looking up the right column/row combos)
  # ------------------------------------------------------------------------------------------
  # *** Stochastic interventions ***
  #     Need to do either MC integration (sample from g.star then predict Q)
  #     or direct weighted sum of predicted Q's with weights given by g.star (on A and N)
  # ------------------------------------------------------------------------------------------
  # 1. At each t iteration:
  #     Fitting:
  #        Outcome is always Q.kplus1 (either by regimen or not).
  #        Need to remove all censored at t when fitting (add new expr to subset defn or pass second subset expr?).
  #        One way is to assign Q.kplus1[t] <- NA to all observations currently censored (at t)
  #        Then use evalsubst() with outvar for fit() and without outvar for predict().
  #     Prediction:
  #        Add observation which were censored at t (C[t]==1) to the prediction set, using the second subset (or outvar)
  #        For each Q.m (Q fit), rename columns to replace current A[t] with A^*[t], (possibly for each regimen), put A[t] back afterwards.
  #        Save ProbA1 from PredictP1() for all obs used in prediction in rows of column Q.kplus1 (either at row Q.kplus1[t] or at Q.kplus1[t-1])
  # 2. At next iteration t-1:
  #     Fitting:
  #        If previous prediction was saved at Q.kplus1[t-1], then all outcomes are already set to their needed values
  #        If not, the outcomes at Q.kplus1[t-1] need to be updated with values from Q.kplus1[t] for all IDs who were used in prediction at t.
  # ------------------------------------------------------------------------------------------

  # Run all Q-learning regressions (one for each subsets defined above, predictions of the last regression form the outcomes for the next:
  Qlearn.fit$fit(data = OData)

  browser()
  OData$dat.sVar[1:50,]
  Qlearn.fit

  # 1a. Grab the mean prediction from the very last regression (over all n observations);
  # this is the G-comp estimate of survival for the initial t value (t used in the first Q-reg)
    lastQ_inx <- Qreg_idx[1] # the index for the last Q-fit
    reslastQP1 <- Qlearn.fit$predictRegK(lastQ_inx, OData$nuniqueIDs)
    print(length(reslastQP1))
    print(mean(reslastQP1))
    # [1] 0.143 under g.0
    # [1] 0.176548 for abar = 000000

  # 1b. Grab the right model (QlearnModel) and pull it directly:
    lastQ.fit <- Qlearn.fit$getPsAsW.models()[[lastQ_inx]]$getPsAsW.models()[[1]]
    lastQ.fit
    # for all observations in long format:
    length(lastQ.fit$getprobA1)
    head(lastQ.fit$getprobA1)
    # for all unique IDs in the data:
    length(lastQ.fit$predictAeqa())
    head(lastQ.fit$predictAeqa())
    mean(lastQ.fit$predictAeqa())

  # 1c. Grab it directly from the data, using the appropriate strata-subsetting expression
    subset_vars <- lastQ.fit$subset_vars
    subset_exprs <- lastQ.fit$subset_exprs
    subset_idx <- OData$evalsubst(subset_vars = subset_vars, subset_exprs = subset_exprs)
    mean(OData$dat.sVar[subset_idx, ][["Q.kplus1"]])

  OData$Qlearn.fit <- Qlearn.fit

  return(OData)
}