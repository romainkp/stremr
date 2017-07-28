## ----------------------------------------------------------------------------------
## Class for defining, fitting and predicting for a single regression model E(Y|X) (univariate outcome).
## R6 class for fitting and making predictions for a single outcome regression model.
## This R6 class can request, store and manage the design matrix Xmat, as well as the outcome Y.
## ----------------------------------------------------------------------------------
BinaryOutcomeModel  <- R6Class(classname = "BinaryOutcomeModel",
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    classify = TRUE,
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    cont.sVar.flag = logical(),
    bw.j = numeric(),
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

        opt_params <- model_contrl[["opt_params"]]
        model_contrl[["opt_params"]] <- NULL

        if (!("estimator" %in% names(opt_params))) opt_params[["estimator"]] <- model_contrl[["estimator"]]
        if (!("family" %in% names(opt_params))) opt_params[["family"]] <- "quasibinomial"
        if (!("distribution" %in% names(opt_params))) opt_params[["distribution"]] <- "bernoulli"

        self$models <- do.call(gridisl::defModel, opt_params)

      }

      self$model_contrl <- model_contrl

      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar

      assert_that(is.character(reg$predvars))
      self$predvars <- reg$predvars

      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)

      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- reg$nbins

      if (is.null(reg$subset_vars)) {self$subset_vars <- TRUE}
      assert_that(is.logical(self$subset_vars) || is.character(self$subset_vars)) # is.call(self$subset_vars) ||

      # if (gvars$verbose) {
      #   print("New instance of " %+% class(self)[1] %+% " :"); print(self$show())
      # }

      # Get the bin width (interval length) for the current bin name self$getoutvarnm (for discretized continuous sA only):
      self$cont.sVar.flag <- self$getoutvarnm %in% names(reg$intrvls.width)
      if (self$cont.sVar.flag) {
        intrvl.idx <- which(names(reg$intrvls.width) %in% self$getoutvarnm)
        if (length(intrvl.idx) > 1) stop("non-unique names for intrvls.width in RegressionClass")
        self$bw.j <- reg$intrvls.width[intrvl.idx]
      } else {
        self$bw.j <- 1L
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
      private$model.fit <- fit_single_regression(data, nodes, self$models, self$model_contrl, self$predvars, self$outvar, self$subset_idx)

      self$is.fitted <- TRUE
      if (predict) self$predictAeqa(..., indA = data$get.outvar(self$subset_idx, self$getoutvarnm))

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # **********************************************************************
      self$wipe.alldat
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$model.fit to generate predictions for data:
    predict = function(newdata, holdout = FALSE, ...) {
      assert_that(self$is.fitted)

      if (missing(newdata) && is.null(private$probA1)) {
        private$probA1 <- gridisl::predict_SL(modelfit = private$model.fit,
                                              add_subject_data = FALSE,
                                              subset_idx = self$subset_idx,
                                              # use_best_retrained_model = TRUE,
                                              holdout = holdout,
                                              verbose = gvars$verbose)

      } else {

        self$n <- newdata$nobs
        self$define.subset.idx(newdata)
        private$probA1 <- gridisl::predict_SL(modelfit = private$model.fit,
                                              newdata = newdata,
                                              add_subject_data = FALSE,
                                              subset_idx = self$subset_idx,
                                              # use_best_retrained_model = TRUE,
                                              holdout = holdout,
                                              verbose = gvars$verbose)

      }
      return(invisible(self))
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, bw.j.sA_diff, indA, holdout = FALSE, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      }

      self$predict(newdata, holdout)

      if (missing(newdata) & missing(indA)) {
        indA <- self$getoutvarval
      } else if (missing(indA)) {
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      }

      assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      probA1 <- private$probA1 # [self$getsubset]

      # check that predictions P(A=1 | dmat) exist for all obs (not NA)
      if (any(is.na(probA1))) {
        stop("some of the modeling predictions resulted in NAs, which indicates an error in a prediction routine")
      }

      # Discrete version for joint density:
      # likelihood_1 <- as.vector(probA1^(indA) * (1 - probA1)^(1L - indA))
      likelihood_2 <- probA1[[1]]^(indA) * (1 - probA1[[1]])^(1L - indA)
      # likelihood_3 <- probA1[, (names(probA1)) := .SD^(indA) * (1 - .SD)^(1L - indA)]

      # probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      probAeqa[self$getsubset] <- likelihood_2

      # Continuous version for the joint density:
      # probAeqa[self$getsubset] <- (probA1^indA) * exp(-probA1)^(1 - indA)
      # Alternative intergrating the last hazard chunk up to x:
      # difference of sA value and its left most bin cutoff: x - b_{j-1}
      if (!missing(bw.j.sA_diff)) {
        # + integrating the constant hazard all the way up to value of each sa:
        # probAeqa[self$getsubset] <- probAeqa[self$getsubset] * (1 - bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
        # cont. version of above:
        probAeqa[self$getsubset] <- probAeqa[self$getsubset] * exp(-bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1[[1]])^(indA)
      }
      private$probAeqa <- probAeqa

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      return(probAeqa)
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

    # take fitted BinaryOutcomeModel class object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("BinaryOutcomeModel" %in% class(bin.out.model))
      private$model.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # take BinaryOutcomeModel class object that contains the predictions for P(A=1|sW) and save these predictions to self$
    copy.predict = function(bin.out.model) {
      assert_that("BinaryOutcomeModel" %in% class(bin.out.model))
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
        return("P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=", ") %+% ")" %+% ";\\ Stratify: " %+% self$subset_exprs %+% ";\\ N: " %+% self$n_obs_fit)
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs, N = self$n_obs_fit))
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