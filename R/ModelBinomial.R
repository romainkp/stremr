## ----------------------------------------------------------------------------------
## Class for defining, fitting and predicting for a single regression model E(Y|X) (Binomial outcome).
## R6 class for fitting and making predictions for a single outcome regression model.
## This R6 class can request, store and manage the design matrix Xmat, as well as the outcome Y.
## ----------------------------------------------------------------------------------
ModelBinomial  <- R6Class(classname = "ModelBinomial",
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    outvar = character(),   # outcome name(s)
    outvar.class = character(),
    outcome_type = character(),
    predvars = character(), # names of predictor vars
    is.fitted = FALSE,
    model_contrl = list(),
    models = list(),
    n = NA_integer_,         # total number of rows in the input data
    n_obs_fit = NA_integer_, # total number of observations used for fitting the model
    nbins = integer(),
    subset_vars = NULL,      # THE VAR NAMES WHICH WILL BE TESTED FOR MISSINGNESS AND WILL DEFINE SUBSETTING
    subset_exprs = NULL,     # THE LOGICAL EXPRESSION (ONE) TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,       # Logical vector of length n (TRUE = include the obs)
    ReplMisVal0 = logical(),

    initialize = function(reg, ...) {
      model_contrl <- reg$model_contrl

      if (!is.null(model_contrl[["models"]])) {
        self$models <- model_contrl[["models"]]
        model_contrl[["models"]] <- NULL
      } else {
        stop("for binomial exposures 'model' must be always always specified and it must be an sl3 learner object")
      }

      self$model_contrl <- model_contrl

      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar
      self$outvar.class <- reg$outvar.class
      self$outcome_type <- "quasibinomial"

      assert_that(is.character(reg$predvars))
      self$predvars <- reg$predvars

      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)

      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- reg$nbins

      if (is.null(reg$subset_vars)) {self$subset_vars <- TRUE}
      assert_that(is.logical(self$subset_vars) || is.character(self$subset_vars)) # is.call(self$subset_vars) ||

      if (gvars$verbose) {
        print("New 'ModelBinomial' regression defined:"); print(self$show())
      }

      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and private$probAeqa
    fit = function(overwrite = FALSE, data, predict = FALSE, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      if (gvars$verbose) print("fitting the model: " %+% self$show())
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      self$define.subset.idx(data)
      nodes <- data$nodes

      self$n_obs_fit <- length(self$subset_idx)

      ## only perform fit if the strata (current training data) is not empty
      if (self$n_obs_fit > 0) {
        private$model.fit <- fit_single_regression(data = data, 
                                                   nodes = nodes, 
                                                   models = self$models, 
                                                   model_contrl = self$model_contrl, 
                                                   predvars = self$predvars, 
                                                   outvar = self$outvar, 
                                                   subset_idx = self$subset_idx,
                                                   outcome_type = self$outcome_type)        
      }

      self$is.fitted <- TRUE
      if (predict) try(self$predictAeqa(..., indA = data$get.outvar(self$subset_idx, self$getoutvarnm)), silent = TRUE)

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # **********************************************************************
      self$wipe.alldat
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$model.fit to generate predictions for data:
    # if prediction data has 0 observations, return an empty prediction vector
    predict = function(newdata, holdout = FALSE, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata) && !is.null(private$probA1)) {
        return(invisible(self))
      }
      model.fit <- private$model.fit
      ## change to accomodate case when subset is 0 length
      if (missing(newdata)) {
        if (length(self$subset_idx) == 0L) { ## to accomodate case when subset is 0 length
          private$probA1 = double()
        } else if (is(model.fit, "PredictionStack")) {
          private$probA1 <- gridisl::predict_SL(modelfit = model.fit, add_subject_data = FALSE, subset_idx = self$subset_idx, holdout = holdout, verbose = gvars$verbose)
        } else if (is(model.fit, "Lrnr_base")) {
          private$probA1 <- model.fit$predict()
        } else {
          stop("model fit object is of unrecognized class (private$model.fit)")
        }
      } else {
        self$n <- newdata$nobs
        self$define.subset.idx(newdata)
        if (length(self$subset_idx) == 0L) { ## to accomodate case when subset is 0 length
          private$probA1 = double()
        } else if (is(model.fit, "PredictionStack")) {
          private$probA1 <- gridisl::predict_SL(modelfit = model.fit, newdata = newdata, add_subject_data = FALSE, subset_idx = self$subset_idx, holdout = holdout, verbose = gvars$verbose)
        } else if (is(model.fit, "Lrnr_base")) {
          ## todo: allow passing the DataStorage object directly to task, seamlessly
          new_task <- sl3::sl3_Task$new(newdata$dat.sVar[self$subset_idx, ], covariates = self$predvars, outcome = self$outvar)
          private$probA1 <- model.fit$predict(new_task)
        } else {
          stop("model fit object is of unrecognized class (private$model.fit)")
        }

      }
      return(invisible(self))
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, indA, holdout = FALSE, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      }
      self$predict(newdata, holdout)

      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      probA1 <- private$probA1 # [self$getsubset]
      # check that predictions P(A=1 | dmat) exist for all obs (not NA)
      if (any(is.na(probA1))) stop("some of the modeling predictions resulted in NAs, which indicates an error in a prediction routine")
      ## 2) regular classification / binary regression problem, turn into likelihood:
      if (missing(newdata) & missing(indA)) {
        indA <- self$getoutvarval
      } else if (missing(indA)) {
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      }
      assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
      if (is.list(probA1) || is.data.table(probA1) || is.data.frame(probA1)) {
        probA1 <- probA1[[1]]
      }
      likelihood <- probA1^(indA) * (1 - probA1)^(1L - indA)

      probAeqa[self$getsubset] <- likelihood
      private$probAeqa <- probAeqa

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      return(probAeqa)
    },

    predictgstar = function(newdata, intervened_NODE_all, intervened_type_all, useonly_t_NODE_all, ...) { # P(A^*[i]=A^*|W=w) - calculating the likelihood for g^*(A)
      assert_that(self$is.fitted)
      self$n <- newdata$nobs
      self$define.subset.idx(newdata)

      if (missing(newdata)) {
        stop("newdata must be provided for evaluating gstar")
      }

      gstar_name = intervened_NODE_all[[self$getoutvarnm]]
      intervened_type = intervened_type_all[[self$getoutvarnm]]
      useonly_t = useonly_t_NODE_all[[gstar_name]]
      # print("useonly_t_NODE_all"); print(useonly_t_NODE_all)
      # print(paste0("useonly_t: ", useonly_t))

      # Aobs <-  newdata$get.outvar(TRUE, self$getoutvarnm)
      # Astar <- newdata$get.outvar(TRUE, gstar_name)
      Aobs <-  newdata$get.outvar(self$getsubset, self$getoutvarnm)
      Astar <- newdata$get.outvar(self$getsubset, gstar_name)

      # gstar <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      if (intervened_type %in% "bin") {
        # -------------------------------------------------
        # map gstar for intervention on binary A
        # -------------------------------------------------
        # check that observed exposure is always a vector of integers
        if(!is.integerish(Aobs)) {
          stop("Not possible to intervene on continuous node ('" %+% self$getoutvarnm %+% 
               "') with a static intervention -- please change the intervention type by setting 'intervened_type_TRT' / 'intervened_type_MONITOR' to 'shift' or 'MSM'.")
        }
        gstar <- Astar^(Aobs) * (1 - Astar)^(1L - Aobs)
      } else if (intervened_type %in% "shift") {
        # -------------------------------------------------
        ## map gstar for a delta(W) shift of continuous A
        # -------------------------------------------------
        ## 1. define the new outcome node to evalute g(A-\delta(W)), defined as Anew = 2*Aobs - Astar
        newdata$dat.sVar[, "Anew.delta.shift.star.tmp" := 2*Aobs - Astar]
        ## 2. swap the names of gnode and self$gstar_name
        newdata$swapNodes(current = self$getoutvarnm, target = "Anew.delta.shift.star.tmp")
        ## 3. call predict on original g fit, but using Anew.tmp as the outcome
        gstar <- self$predictAeqa(newdata)
        ## 4. swap the node names back
        newdata$swapNodes(current = "Anew.delta.shift.star.tmp", target = self$getoutvarnm)
        newdata$dat.sVar[, "Anew.delta.shift.star.tmp" := NULL]
      } else if (intervened_type %in% "MSM") {
        # -------------------------------------------------
        ## map gstar for MSMs, which is just const 1
        # -------------------------------------------------
        gstar <- rep.int(1L, length(self$getsubset))
      } else {
        stop("Unrecognized value for 'intervened_type_TRT' or 'intervened_type_MONITOR', these can only be 'bin', 'shift' or 'MSM', please consult the manual.")
      }

      # print(paste0("Evaluating g^* for exposure node: ", 
      #         self$getoutvarnm, " and intervention node: ", gstar_name))
      # print(paste0("Evaluating g^* for intervened_type: ", intervened_type))
      # print(paste0("Length subset: ", length(self$getsubset)))
      # print("...Aobs..."); print(head(Aobs)); 
      # print("...Astar..."); print(head(Astar)); 
      # print("...gstar..."); print(head(gstar)); 
      # print(paste0("Aobs  len: ", length(Aobs), "; sum: ", sum(Aobs))); 
      # print(paste0("Astar len: ", length(Astar), "; sum: ", sum(Astar)));
      # print(paste0("gstar len: ", length(gstar), "; sum: ", sum(gstar)));

      ## Changing the default assignment value for g^*(A) to NA
      gstar_out <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      # gstar_out <- rep.int(NA_real_, self$n) # this will not work when there is more than one strata to fit g
      gstar_out[self$getsubset] <- gstar

      ## only intervene (assign to gstar) on observations in this subset
      subset_idx_intervene <- newdata$evalsubst(subset_exprs = useonly_t)
      if (any(is.na(subset_idx_intervene)))
        stop("the subset index evaluation for the expression '" %+% useonly_t %+% "' resulted in NAs")
      idx_set_to_g0 <- setdiff(self$getsubset, subset_idx_intervene)
      gstar_out[idx_set_to_g0] <- private$probAeqa[idx_set_to_g0]

      self$wipe.alldat # to save RAM space when doing many stacked regressions wipe out all internal data:

      return(gstar_out)
    },

    define.subset.idx = function(data) {
      if (is.logical(self$subset_vars)) {
        subset_idx <- which(self$subset_vars)
      } else if (is.call(self$subset_vars)) {
        stop("calls aren't allowed for self$subset_vars")
      } else if (is.character(self$subset_vars)) {
        subset_idx <- data$evalsubst(subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)
      }
      self$subset_idx <- subset_idx
      return(invisible(self))
    },

    # take fitted ModelBinomial class object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("ModelBinomial" %in% class(bin.out.model))
      private$model.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # take ModelBinomial class object that contains the predictions for P(A=1|sW) and save these predictions to self$
    copy.predict = function(bin.out.model) {
      assert_that("ModelBinomial" %in% class(bin.out.model))
      assert_that(self$is.fitted)
      private$probA1 <- bin.out.model$getprobA1
    },

    # Returns the object that contains the actual model fits (itself)
    get.fits = function() {
      model.fit <- self$getfit
      return(list(model.fit))
    },

    get.model.summaries = function() {
      return(list(self$show(print_format = FALSE)))
    },

    # Output info on the general type of regression being fitted:
    show = function(print_format = TRUE) {
      if (print_format) {
        return("P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=", ") %+% ")" %+% 
              ";\\ outvar.class: " %+% self$outvar.class %+% 
              ";\\ Stratify: " %+% self$subset_exprs %+% 
              ";\\ N: " %+% self$n_obs_fit)
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs, N = self$n_obs_fit, outvar.class = self$outvar.class))
      }
    }
  ),

  active = list(
    wipe.alldat = function() {
      # private$probA1 <- NULL
      # private$probAeqa <- NULL
      self$subset_idx <- NULL
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    getoutvarval = function() { stop("self$getoutvarval is not implemented") }
  ),
  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    probA1 = NULL,    # Predicted probA^s=1 conditional on Xmat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)