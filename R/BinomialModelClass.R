#----------------------------------------------------------------------------------
# Classes for modelling regression models with binary outcome Bin ~ Xmat
#----------------------------------------------------------------------------------

## ---------------------------------------------------------------------
#' R6 class for fitting and making predictions for a single binary outcome regression model P(B | PredVars)
#'
#' This R6 class can request, store and manage the design matrix Xmat, as well as the binary outcome Bin for the
#'  logistic regression P(Bin|Xmat).
#'  Can also be used for converting data in wide format to long when requested,
#'  e.g., when pooling across binary indicators (fitting one pooled logistic regression model for several indicators)
#'  The class has methods that perform queries to data storage R6 class DataStorageClass to get appropriate data columns & row subsets
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{cont.sVar.flag} - Is the original outcome variable continuous?
#' \item{bw.j} - Bin width (interval length) for an outcome that is a bin indicator of a discretized continous outcome.
#' \item{GLMpackage} - Controls which package will be used for performing model fits (\code{glm} or \code{speedglm}).
#' \item{binomialModelObj} - Pointer to an instance of \code{binomialModelObj} class that contains the data.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg)}}{Uses \code{reg} R6 \code{\link{RegressionClass}} object to instantiate a new model for a
#'   logistic regression with binary outcome.}
#'   \item{\code{show()}}{Print information on outcome and predictor names used in this regression model}
#'   \item{\code{fit()}}{...}
#'   \item{\code{copy.fit()}}{...}
#'   \item{\code{predict()}}{...}
#'   \item{\code{copy.predict()}}{...}
#'   \item{\code{predictAeqa()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{getoutvarnm}}{...}
#'   \item{\code{getoutvarval}}{...}
#'   \item{\code{getsubset}}{...}
#'   \item{\code{getprobA1}}{...}
#'   \item{\code{getfit}}{...}
#'   \item{\code{wipe.alldat}}{...}
#' }
#' @importFrom assertthat assert_that is.flag
#' @export
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

    binomialModelObj = NULL, # object of class binomialModelObj that is used in fitting / prediction, never saved (need to be initialized with $new())
    fit.package = character(),
    fit.algorithm = character(),
    model_contrl = list(),

    n = NA_integer_,        # number of rows in the input data
    nbins = integer(),
    subset_vars = NULL,     # THE VAR NAMES WHICH WILL BE TESTED FOR MISSINGNESS AND WILL DEFINE SUBSETTING
    subset_exprs = NULL,     # THE LOGICAL EXPRESSION (ONE) TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)

    ReplMisVal0 = logical(),

    initialize = function(reg, ...) {
      self$model_contrl <- reg$model_contrl

      if ("fit.package" %in% names(self$model_contrl)) {
        self$fit.package <- self$model_contrl[['fit.package']]
        assert_that(is.character(self$fit.package))
      } else {
        self$fit.package <- reg$fit.package[1]
      }
      if (!(self$fit.package %in% allowed.fit.package)) stop("fit.package must be one of: " %+% paste0(allowed.fit.package, collapse=", "))

      if ("fit.algorithm" %in% names(self$model_contrl)) {
        self$fit.algorithm <- self$model_contrl[['fit.algorithm']]
        assert_that(is.character(self$fit.algorithm))
      } else {
        self$fit.algorithm <- reg$fit.algorithm[1]
      }
      if (!(self$fit.algorithm %in% allowed.fit.algorithm)) stop("fit.algorithm must be one of: " %+% paste0(allowed.fit.algorithm, collapse=", "))

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

      # ***************************************************************************
      # Add any additional options passed on to modeling functions as extra args
      # ***************************************************************************
      if (self$fit.package %in% c("h2o", "h2oEnsemble")) {
        self$binomialModelObj <- BinomialH2O$new(fit.algorithm = self$fit.algorithm, fit.package = self$fit.package, ParentModel = self, ...)
      } else {
        self$binomialModelObj <- BinomialGLM$new(fit.algorithm = self$fit.algorithm, fit.package = self$fit.package, ParentModel = self, ...)
      }
      self$binomialModelObj$params <- list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs)

      if (gvars$verbose) {
        print("New instance of " %+% class(self)[1] %+% " :"); print(self$show())
      }
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
      model.fit <- self$binomialModelObj$fit(data, self$outvar, self$predvars, self$subset_idx, ...)

      if (inherits(model.fit, "try-error")) {
        message("running " %+% self$binomialModelObj$fit.class %+% " with h2o has failed, trying to run speedglm as a backup...")
        self$binomialModelObj <- BinomialGLM$new(fit.algorithm = "glm", fit.package = "speedglm", ParentModel = self, ...)
        self$binomialModelObj$params <- list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs)
        model.fit <- self$binomialModelObj$fit(data, self$outvar, self$predvars, self$subset_idx, ...)
      }

      private$model.fit <- model.fit

      self$is.fitted <- TRUE
      if (predict) {
        self$predictAeqa(...)
      }

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # **********************************************************************
      self$wipe.alldat
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$model.fit to generate predictions for data:
    predict = function(newdata, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata) && is.null(private$probA1)) {
        private$probA1 <- self$binomialModelObj$predictP1(subset_idx = self$subset_idx)
      } else {
        self$n <- newdata$nobs
        self$define.subset.idx(newdata)
        private$probA1 <- self$binomialModelObj$predictP1(data = newdata, subset_idx = self$subset_idx)
      }
      return(invisible(self))
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, bw.j.sA_diff, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      }

      self$predict(newdata)

      if (missing(newdata)) {
        indA <- self$getoutvarval
      } else {
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      }

      assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      probA1 <- private$probA1[self$getsubset]

      # check that predictions P(A=1 | dmat) exist for all obs (not NA)
      if (any(is.na(probA1) & !is.nan(probA1))) {
        stop("some of the modeling predictions resulted in NAs, which indicates an error of a prediction routine")
      }
      # assert_that(!any(is.na(probA1)))

      # Discrete version for joint density:
      probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      # continuous version for the joint density:
      # probAeqa[self$getsubset] <- (probA1^indA) * exp(-probA1)^(1 - indA)
      # Alternative intergrating the last hazard chunk up to x:
      # difference of sA value and its left most bin cutoff: x - b_{j-1}
      if (!missing(bw.j.sA_diff)) {
        # + integrating the constant hazard all the way up to value of each sa:
        # probAeqa[self$getsubset] <- probAeqa[self$getsubset] * (1 - bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
        # cont. version of above:
        probAeqa[self$getsubset] <- probAeqa[self$getsubset] * exp(-bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
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
        subset_idx <- self$subset_vars
      } else if (is.call(self$subset_vars)) {
        stop("calls aren't allowed in binomialModelObj$subset_vars")
      } else if (is.character(self$subset_vars)) {
        subset_idx <- data$evalsubst(subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)
      }
      assert_that(is.logical(subset_idx))
      if ((length(subset_idx) < self$n) && (length(subset_idx) > 1L)) {
        if (gvars$verbose) message("subset_idx has smaller length than self$n; repeating subset_idx p times, for p: " %+% data$p)
        subset_idx <- rep.int(subset_idx, data$p)
        if (length(subset_idx) != self$n) stop("binomialModelObj$define.subset.idx: self$n is not equal to nobs*p!")
      }
      assert_that((length(subset_idx) == self$n) || (length(subset_idx) == 1L))
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
      # private$probA1 <- NULL
      # private$probAeqa <- NULL
      self$subset_idx <- NULL
      self$binomialModelObj$emptydata
      self$binomialModelObj$emptyY
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    getoutvarval = function() { self$binomialModelObj$getY }
  ),
  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    probA1 = NULL,    # Predicted probA^s=1 conditional on Xmat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)

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
      self$subset_idx <- rep.int(TRUE, self$n)
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
      probA1 <- private$probA1[self$getsubset]
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