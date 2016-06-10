#----------------------------------------------------------------------------------
# Classes for modelling regression models with binary outcome Bin ~ Xmat
#----------------------------------------------------------------------------------

# S3 methods for getting model fit for BinaryOutcomeModel class object
fit.BinaryOutcomeModel <- function(BinaryOutcomeModel) {
  assert_that(BinaryOutcomeModel$is.fitted)
  BinaryOutcomeModel$getfit
}
# S3 methods for getting model fit summary for BinaryOutcomeModel class object
summary.BinaryOutcomeModel <- function(BinaryOutcomeModel) {
  assert_that(BinaryOutcomeModel$is.fitted)
  fit <- BinaryOutcomeModel$getfit
  append(list(reg = BinaryOutcomeModel$show()), fit)
}

## ---------------------------------------------------------------------
#' R6 class for fitting and making predictions for a single logistic regression with binary outcome B, P(B | PredVars)
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
#' \item{glmfitclass} - Controls which package will be used for performing model fits (\code{glm} or \code{speedglm}).
#' \item{bindat} - Pointer to an instance of \code{BinDat} class that contains the data.
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
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    cont.sVar.flag = logical(),
    bw.j = numeric(),
    glmfitclass = "glmS3", # default glm fit class
    is.fitted = FALSE,
    bindat = NULL, # object of class BinDat that is used in fitting / prediction, never saved (need to be initialized with $new())

    initialize = function(reg, ...) {
      assert_that(is.character(reg$GLMpackage))
      if (reg$GLMpackage %in% "glm") {
        self$glmfitclass <- "glmS3"
      } else if (reg$GLMpackage %in% "speedglm") {
        self$glmfitclass <- "speedglmS3"
      } else if (reg$GLMpackage %in% "h2o") {
        self$glmfitclass <- "h2oglmS3"
      } else {
        stop("reg$GLMpackage is of unrecognized type")
      }

      self$outvar <- reg$outvar
      self$predvars <- reg$predvars

      self$bindat <- BinDat$new(reg = reg, ...)
      class(self$bindat) <- c(class(self$bindat), self$glmfitclass)
      if (gvars$verbose) {
        print("New instance of BinaryOutcomeModel:"); print(self$show())
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

    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      if (gvars$verbose) print("fitting the model: " %+% self$show())
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked
      self$bindat$newdata(newdata = data, ...) # populate bindat with X_mat & Y_vals
      private$m.fit <- logisfit(BinDatObject = self$bindat) # private$m.fit <- data_obj$logisfit or private$m.fit <- data_obj$logisfit()
      # alternative 2 is to apply data_obj method / method that fits the model
      self$is.fitted <- TRUE
      self$wipe.alldat
      invisible(self)
    },

    # take fitted BinaryOutcomeModel class object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("BinaryOutcomeModel" %in% class(bin.out.model))
      private$m.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$m.fit to generate predictions for newdata:
    predict = function(newdata, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata)) {
        stop("must provide newdata for BinaryOutcomeModel$predict()")
      }
      # re-populate bindat with new X_mat:
      self$bindat$newdata(newdata = newdata, getoutvar = FALSE, ...)
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        stop("BinaryOutcomeModel$predict is not applicable to pooled regression, call BinaryOutcomeModel$predictAeqa instead")
      } else {
        private$probA1 <- logispredict(m.fit = private$m.fit, BinDatObject = self$bindat)
      }
      self$bindat$emptydata  # Xmat in bindat is no longer needed, but subset, outvar & probA1 may be needed for private$probA1
      invisible(self)
    },

    # take BinaryOutcomeModel class object that contains the predictions for P(A=1|sW) and save these predictions to self$
    copy.predict = function(bin.out.model) {
      assert_that("BinaryOutcomeModel" %in% class(bin.out.model))
      assert_that(self$is.fitted)
      private$probA1 <- bin.out.model$getprobA1
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, bw.j.sA_diff) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      assert_that(!missing(newdata))
      self$bindat$newdata(newdata = newdata, getoutvar = TRUE) # populate bindat with new design matrix covars X_mat
      assert_that(is.logical(self$getsubset))
      n <- newdata$nobs
      # obtain predictions (likelihood) for response on fitted data (from long pooled regression):
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        probAeqa <- logispredict.long(m.fit = private$m.fit, BinDatObject = self$bindat) # overwrite probA1 with new predictions:
      } else {
        # get predictions for P(A[j]=1|W=newdata) from newdata:
        probA1 <- logispredict(m.fit = private$m.fit, BinDatObject = self$bindat)
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
        assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
        probAeqa <- rep.int(1L, n) # for missing, the likelihood is always set to P(A = a) = 1.
        assert_that(!any(is.na(probA1[self$getsubset]))) # check that predictions P(A=1 | dmat) exist for all obs.
        probA1 <- probA1[self$getsubset]
        # discrete version for the joint density:
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
      }
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      return(probAeqa)
    },

    sampleA = function(newdata, bw.j.sA_diff) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      assert_that(!missing(newdata))
      # Don't want to subset by the outvar, since binarized mat for cat outcome is not re-created when just sampling
      # But need to reset it back when done
      temp_subset_vars <- self$bindat$subset_vars
      self$bindat$subset_vars <- self$bindat$subset_vars[!self$bindat$subset_vars %in% self$bindat$outvar]
      self$bindat$newdata(newdata = newdata, getoutvar = FALSE) # populate bindat with new design matrix covars X_mat

      assert_that(is.logical(self$getsubset))
      n <- newdata$nobs
      # obtain predictions (likelihood) for response on fitted data (from long pooled regression):
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        stop("not implemented")
      } else {
        # get probability P(sA[j]=1|sW=newdata) from newdata, then sample from rbinom
        probA1 <- logispredict(m.fit = private$m.fit, BinDatObject = self$bindat)
        sampleA <- rep.int(0L, n)
        sampleA[self$getsubset] <- rbinom(n = n, size = 1, prob = probA1)
      }
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      self$bindat$subset_vars <- temp_subset_vars
      # **********************************************************************
      return(sampleA)
    },

    get.fits = function(format_table = FALSE) {
      coef_out <- private$m.fit$coef
      nobs <- private$m.fit$nobs
      if (format_table) {
        coef_out <- data.frame(Coef = coef_out)
        coef_out <- cbind(names(private$m.fit$coef), coef_out)
        rownames(coef_out) <- NULL
        colnames(coef_out) <- c("Terms", "Coefficients")
      }
      return(list(c(self$show(print_format = FALSE), list(nobs = nobs, coef = coef_out))))
    },
    show = function(print_format = TRUE) {self$bindat$show(print_format)}
  ),
  active = list(
    wipe.alldat = function() {
      private$probA1 <- NULL
      private$probAeqa <- NULL
      self$bindat$emptydata
      self$bindat$emptyY
      self$bindat$emptySubset_idx
      self$bindat$emptyN
      return(self)
    },
    getfit = function() { private$m.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$bindat$subset_idx },
    getoutvarval = function() { self$bindat$getY },
    getoutvarnm = function() { self$bindat$outvar }
  ),
  private = list(
    m.fit = list(),   # the model fit (either coefficients or the model fit object)
    probA1 = NULL,    # Predicted probA^s=1 conditional on X_mat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on X_mat
  )
)
