
## ----------------------------------------------------------------------------------
## Call \code{gridisl} and fit a single regression model.
## All model fitting in stremr is performed via this function.
## ----------------------------------------------------------------------------------
fit_single_regression <- function(data, nodes, models, model_contrl, predvars, outvar, subset_idx, ...) {

  if (is.null(model_contrl[["fit_method"]]))
    stop("'fit_method' must be specified")

  method <- model_contrl[["fit_method"]]
  fold_column <- model_contrl[["fold_column"]]

  if ((method %in% c("cv", "origamiSL")) && is.null(fold_column) && is.null(data$fold_column)) {

    stop(
"Must manually define the folds and specify the 'fold_column' to be able to perform cross-validation (method = 'cv').
The fold column can be defined by either:
a) Calling define_CVfolds(data, ...), where data is the object returned by importData(); or
b) Setting the 'fold_column' option to the name of the column that contains the fold IDs (stremrOptions('fold_column', 'name_of_the_fold_column')); or
c) Passing the name of the existing fold column as the argument 'fold_column' of the calling function")

  ## Manually provided column name with fold IDs, perform some check and save the fold information
  } else if ((method %in% c("cv", "origamiSL")) && !is.null(fold_column)) {

    if (!is.character(fold_column)) stop("argument 'fold_column' must be a string")
    data$define_CVfolds(fold_column = fold_column)

  ## Use the existing fold ID column (previously defined by calling define_CVfolds())
  } else if ((method %in% c("cv", "origamiSL")) && is.null(fold_column)) fold_column <- data$fold_column

  model.fit <- try({
    gridisl::fit(models,
                method = method,
                ID = nodes$IDnode,
                t_name = nodes$tnode,
                x = predvars,
                y = outvar,
                data = data,
                fold_column = fold_column,
                subset_idx = subset_idx, ...)
  })

  if (inherits(model.fit, "try-error")) {
    message("...trying to run speedglm as a backup...")
    method <- "none"
    # model_contrl[["fit.package"]] <- "speedglm"
    # model_contrl[["fit.algorithm"]] <- "glm"
    glm_model <- models[1]
    glm_model[[1]][["fit.package"]] <- "speedglm"
    glm_model[[1]][["fit.algorithm"]] <- "glm"
    class(glm_model) <- c(class(glm_model), "ModelStack")
    # glm_model <- gridisl::defModel(estimator = "speedglm__glm", family = family, distribution = distribution)

    model.fit <- gridisl::fit(glm_model,
                             method = method,
                             ID = nodes$IDnode,
                             t_name = nodes$tnode,
                             x = predvars,
                             y = outvar,
                             data = data,
                             verbose = gvars$verbose,
                             fold_column = fold_column,
                             subset_idx = subset_idx)

  }

  return(model.fit)
}

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
      # self$binomialModelObj$emptydata
      # self$binomialModelObj$emptyY
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    # getoutvarval = function() { self$binomialModelObj$getY }
    getoutvarval = function() { stop("self$getoutvarval is not implemented") }
  ),
  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    probA1 = NULL,    # Predicted probA^s=1 conditional on Xmat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)

## ----------------------------------------------------------------------------------
## A trivial class for dealing with deterministic outcome modeling
## ----------------------------------------------------------------------------------
DeterministicBinaryOutcomeModel  <- R6Class(classname = "DeterministicBinaryOutcomeModel",
  inherit = BinaryOutcomeModel,
  cloneable = TRUE,
  portable = TRUE,
  class = TRUE,
  public = list(
    gstar.Name = character(),
    is.fitted = TRUE,

    initialize = function(reg, ...) {
      self$model_contrl <- reg$model_contrl
      self$gstar.Name <- reg$model_contrl[["gstar.Name"]]
      assert_that(!is.null(self$gstar.Name))
      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)
      self$ReplMisVal0 <- reg$ReplMisVal0
      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and private$probAeqa
    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$define.subset.idx(data)
      private$probA1 <- data$get.outvar(TRUE, self$gstar.Name)
      # private$.isNA.probA1 <- is.na(private$probA1)
      # self$subset_idx <- rep.int(TRUE, self$n)
      self$subset_idx <- seq_len(self$n)
      private$.outvar <- data$get.outvar(TRUE, self$getoutvarnm) # Always a vector of 0/1
      # private$.isNA.outvar <- is.na(private$.outvar)
      self$is.fitted <- TRUE
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # self$wipe.alldat
      # **********************************************************************
      invisible(self)
    },

    # get the fixed (known) the gstar P(A^*(t) = 1|W, bar{L(t)});
    # should be already saved earlier in private$probA1, so there is nothing to do here
    predict = function(newdata, ...) {
      assert_that(self$is.fitted)
      return(invisible(self))
    },

    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      if (missing(newdata)) {
        indA <- self$getoutvarval
      } else {
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      }
      assert_that(is.integerish(indA)) # check that observed exposure is always a vector of integers
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      # probA1 <- private$probA1[self$getsubset]
      probA1 <- private$probA1
      probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      self$wipe.alldat # to save RAM space when doing many stacked regressions wipe out all internal data:
      return(probAeqa)
    },

    # Output info on the general type of regression being fitted:
    show = function(print_format = TRUE) {
      if (print_format) {
        return("P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=", ") %+% ")" %+% ";\\ Stratify: " %+% self$subset_exprs)
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs))
      }
    }
  ),

  active = list(
    wipe.alldat = function() {
      private$probA1 <- NULL
      private$probAeqa <- NULL
      private$.outvar <- NULL
      self$subset_idx <- NULL
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    getoutvarval = function() { private$.outvar }
  ),

  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    .outvar = NULL,
    # .isNA.outvar = NULL,
    probA1 = NULL,      # Predicted probA^s=1 conditional on Xmat
    # .isNA.probA1 = NULL,
    probAeqa = NULL     # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)

## ----------------------------------------------------------------------------------
## A trivial class for dealing with NULL outcome modeling (when MONITOR and / or CENS aren't specified)
## This class does nothing but simply returns a vector of (1,1,1,...) when predict methods are called.
## ----------------------------------------------------------------------------------
NULLOutcomeModel  <- R6Class(classname = "NULLOutcomeModel",
  inherit = DeterministicBinaryOutcomeModel,
  cloneable = TRUE,
  portable = TRUE,
  class = TRUE,
  public = list(
    # gstar.Name = character(),
    is.fitted = TRUE,
    initialize = function(reg, ...) {
      self$model_contrl <- reg$model_contrl
      # self$gstar.Name <- reg$model_contrl[["gstar.Name"]]
      # assert_that(!is.null(self$gstar.Name))
      assert_that(is.null(reg$outvar) || reg$outvar == "NULL")
      self$outvar <- reg$outvar
      # self$predvars <- reg$predvars
      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)
      self$ReplMisVal0 <- reg$ReplMisVal0
      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and private$probAeqa
    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      # self$define.subset.idx(data)
      # private$probA1 <- rep.int(1L, self$n)
      # private$.isNA.probA1 <- is.na(private$probA1)
      # self$subset_idx <- rep.int(TRUE, self$n)
      # self$subset_idx <- seq_len(self$n)
      # private$.outvar <- data$get.outvar(TRUE, self$getoutvarnm) # Always a vector of 0/1
      # private$.isNA.outvar <- is.na(private$.outvar)
      # self$is.fitted <- TRUE
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # self$wipe.alldat
      # **********************************************************************
      invisible(self)
    },

    # get the fixed (known) the gstar P(A^*(t) = 1|W, bar{L(t)});
    # should be already saved earlier in private$probA1, so there is nothing to do here
    predict = function(newdata, ...) {
      # assert_that(self$is.fitted)
      return(invisible(self))
    },

    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      # assert_that(self$is.fitted)
      # if (missing(newdata)) {
      #   indA <- self$getoutvarval
      # } else {
      #   indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      # }
      # assert_that(is.integerish(indA)) # check that observed exposure is always a vector of integers
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      # probA1 <- private$probA1[self$getsubset]
      # probA1 <- private$probA1
      # probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      self$wipe.alldat # to save RAM space when doing many stacked regressions wipe out all internal data:
      return(probAeqa)
    },

    # Output info on the general type of regression being fitted:
    show = function(print_format = TRUE) {
      if (print_format) {
        return("P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=", ") %+% ")" %+% ";\\ Stratify: " %+% self$subset_exprs)
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs))
      }
    }
  ),

  active = list(
    wipe.alldat = function() {
      private$probA1 <- NULL
      private$probAeqa <- NULL
      private$.outvar <- NULL
      self$subset_idx <- NULL
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    getoutvarval = function() { private$.outvar }
  ),

  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    .outvar = NULL,
    # .isNA.outvar = NULL,
    probA1 = NULL,      # Predicted probA^s=1 conditional on Xmat
    # .isNA.probA1 = NULL,
    probAeqa = NULL     # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)