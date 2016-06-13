#' @import data.table
NULL

# Generic for fitting the logistic (binomial family) GLM model
glmfit <- function(fit, ...) UseMethod("glmfit")

# S3 method for glm binomial family fit, takes BinDat data object:
# glmfit.glmS3 <- function(BinDatObject, ...) {
glmfit.glm <- function(fit, Xmat, Yvals, ...) {
  if (gvars$verbose) print("calling glm.fit...")
  ctrl <- glm.control(trace = FALSE)
  SuppressGivenWarnings({
    model.fit <- stats::glm.fit(x = Xmat,
                                y = Yvals,
                                family = binomial() ,
                                control = ctrl)
  }, GetWarningsToSuppress())

  fit$coef <- model.fit$coef; fit$fitfunname <- "glm"; fit$linkfun <- "logit_linkinv"; fit$nobs <- nrow(Xmat)

  if (gvars$verbose) {
    print("glm fits:")
    print(fit$coef)
  }
  class(fit)[2] <- "glmfit"

  return(fit)
}

# S3 method for speedglm binomial family fit, takes BinDat data object:
glmfit.speedglm <- function(fit, Xmat, Yvals, ...) {
  if (gvars$verbose) print("calling speedglm.wfit...")
  model.fit <- try(speedglm::speedglm.wfit(X = Xmat,
                                           y = Yvals,
                                           family = binomial(),
                                           trace = FALSE),
                      silent = TRUE)
  if (inherits(model.fit, "try-error")) { # if failed, fall back on stats::glm
    message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", model.fit)
    return(glmfit.glm(fit, Xmat, Yvals, ...))
  }

  fit$coef <- model.fit$coef; fit$fitfunname <- "speedglm"; fit$linkfun <- "logit_linkinv"; fit$nobs <- nrow(Xmat)

  if (gvars$verbose) {
    print("speedglm fits:")
    print(fit$coef)
  }
  class(fit)[2] <- "glmfit"

  return(fit)
}

# S3 method for h2o binomial family fit, takes BinDat data object:
glmfit.h2oglm <- function(fit, DataStorageObject, outvar, predvars, subset_idx, ...) {
  if (length(predvars) == 0L) {
    return(glmfit.speedglm(fit, ...))
  }
  rows_subset <- which(subset_idx)
  subsetH2Oframe <- DataStorageObject$H2O.dat.sVar[rows_subset,]
  outfactors <- as.vector(h2o.unique(subsetH2Oframe[, outvar]))
  NAfactors <- any(is.na(outfactors)) # means that conversion H2O.FRAME produced errors, since there are no NAs in the original data

  if (length(outfactors) < 2L | NAfactors) {
    return(glmfit.speedglm(fit, ...))
  } else if (length(outfactors) > 2L) {
    stop("Attempting to run logistic regression for outcome with more than 2 categories")
  }

  if (gvars$verbose) print("calling h2o.glm...")

  model.fit <- try(h2o::h2o.glm(x = predvars,
                                y = outvar,
                                intercept = TRUE,
                                training_frame = subsetH2Oframe,
                                family = "binomial",
                                standardize = TRUE,
                                solver = c("L_BFGS"), # solver = c("IRLSM"),
                                # remove_collinear_columns = TRUE,
                                max_iterations = 50,
                                lambda = 0L),
                   silent = TRUE)
    # print("h2o.glm fit"); print(model.fit)
    # print(model.fit@parameters)
    # print(model.fit@allparameters)
    # str(model.fit@model)
  if (inherits(model.fit, "try-error")) { # if failed, fall back on stats::glm
    message("h2o::h2o.glm failed, falling back on speedglm; ", model.fit)
    return(glmfit.speedglm(fit, ...))
  }

  # assign the fitted coefficients in correct order (same as predictor order in predvars)
  out_coef <- vector(mode = "numeric", length = length(predvars)+1)
  out_coef[] <- NA
  names(out_coef) <- c("Intercept", predvars)
  out_coef[names(model.fit@model$coefficients)] <- model.fit@model$coefficients

  fit$coef <- out_coef; fit$fitfunname <- "h2o.glm"; fit$linkfun <- "logit_linkinv"; fit$nobs <- length(subset_idx);
  fit$H2O.model.object <- model.fit

  if (gvars$verbose) {
    print("h2oglm fits:")
    print(fit$coef)
  }

  # class(fit)[2] <- "glmfit"
  class(fit)[2] <- "h2ofit"
  return(fit)
}

# ----------------------------------------------------------------
# Generic for predicting P(A=1|...) from either logistic (binomial family) GLM model or H2O model
# ----------------------------------------------------------------
predictP1 <- function(m.fit, ...) UseMethod("predictP1")

# ----------------------------------------------------------------
# Generic prediction fun for logistic regression coefs, predicts P(A = 1 | newXmat)
# predictP1.glmfit <- function(m.fit, Xmat, subset_idx, n, ...) {
# ----------------------------------------------------------------
predictP1.glmfit <- function(m.fit, ParentObject, DataStorageObject, subset_idx, n, ...) {
# predictP1.glm <- function(m.fit, ParentObject, subset_idx, n, ...) {
  ParentObject$setdata(DataStorageObject, subset_idx = subset_idx, getoutvar = FALSE, getXmat = TRUE)
  Xmat <- ParentObject$getXmat
  assert_that(!is.null(Xmat)); assert_that(!is.null(subset_idx))
  # Set to default missing value for A[i] degenerate/degerministic/misval:
  # Alternative, set to default replacement val: pAout <- rep.int(gvars$misXreplace, newBinDatObject$n)
  pAout <- rep.int(gvars$misval, n)
  if (sum(subset_idx > 0)) {
    eta <- Xmat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
    pAout[subset_idx] <- match.fun(FUN = m.fit$linkfun)(eta)
  }
  return(pAout)
}
predictP1.speedglm <- function(m.fit, ...) predictP1.glm(m.fit, ...)

# Two ways to perform prediction for GLM with H2O, BOTH SHOULD WORK THE SAME
predictP1.h2oglm <- function(m.fit, ...) predictP1.glm(m.fit, ...)
# predictP1.h2oglm <- function(m.fit, ...) predictP1.h2o(m.fit, ...)

# ----------------------------------------------------------------
# SAME PREDICTION FUNCTION SHOULD APPLY TO ALL H2O FITTED MODELS
# ----------------------------------------------------------------
predictP1.h2ofit <- function(m.fit, ParentObject, DataStorageObject, subset_idx, n, ...) {
  assert_that(!is.null(subset_idx))
  rows_subset <- which(subset_idx)
  subsetH2Oframe <- DataStorageObject$H2O.dat.sVar[rows_subset,]

  ParentObject$setdata(DataStorageObject, subset_idx = subset_idx, getoutvar = FALSE, getXmat = FALSE)
  pAout <- rep.int(gvars$misval, n)
  if (sum(subset_idx > 0)) {
    pAout[subset_idx] <- as.vector(h2o::h2o.predict(m.fit$H2O.model.object, newdata = subsetH2Oframe)[,"p1"])
  }
  return(pAout)
}


## ---------------------------------------------------------------------
#' R6 class for storing the design matrix and the binary outcome for a single GLM (logistic) regression
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
#' \item{ID} - Vector of observation IDs, \code{1:n}, used for pooling.
#' \item{outvar} - Outcome name.
#' \item{predvars} - Predictor names.
#' \item{subset_vars} - Defines the subset which would be used for fitting this model (logical, expression or indices).
#' \item{subset_idx} - Subset \code{subset_vars} converted to logical vector.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg)}}{Uses \code{reg} R6 \code{\link{RegressionClass}} object to instantiate a new storage container for a
#'   design matrix and binary outcome.}
#'   \item{\code{setdata()}}{...
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{emptydata}}{...}
#'   \item{\code{emptyY}}{...}
#'   \item{\code{emptySubset_idx}}{...}
#'   \item{\code{getXmat}}{...}
#'   \item{\code{getY}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.string is.flag
#' @export
# BinDat <- R6Class(classname = "BinDat",
BinomialGLM <- R6Class(classname = "BinomialGLM",
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    ParentModel = NULL,

    # TO DO: THIS WILL CONTAIN ADDITIONAL USER-SPEC'ED CONTROLS/ARGS PASSED ON TO speedglm/glm
    model.controls = NULL,
    fit.class = c("glm", "speedglm", "h2oglm"),
    model.fit = list(coef = NA, fitfunname = NA, nobs = NA),
    # fit.algorithm <- "glm"

    initialize = function(fit.algorithm, fit.package, ParentModel, ...) {
      self$ParentModel <- ParentModel

      if (!("glm" %in% fit.algorithm)) warning("over-riding fit.algorithm option with 'glm', since fit.package was set to 'speedglm' or 'glm' or 'h2oglm'")
      # self$fit.algorithm <- "glm"
      assert_that(any(c("glm", "speedglm", "h2oglm") %in% fit.package))
      self$fit.class <- fit.package
      class(self$model.fit) <- c(class(self$model.fit), self$fit.class)

      invisible(self)
    },

    fit = function(data, outvar, predvars, subset_idx, ...) {
      self$setdata(data, subset_idx = subset_idx, getXmat = TRUE, ...)
      # Xmat has 0 rows: return NA's and avoid throwing exception:
      if (sum(subset_idx) == 0L) {
        self$model.fit$coef = rep.int(NA_real_, ncol(private$Xmat))
      } else {
        self$model.fit <- glmfit(self$model.fit,
                                 Xmat = private$Xmat,
                                 Yvals = private$Yvals,
                                 DataStorageObject = data,
                                 outvar = outvar,
                                 predvars = predvars,
                                 subset_idx = subset_idx, ...)
      }
      return(self$model.fit)
    },

    predictP1 = function(data, subset_idx) {
      # self$setdata(data, subset_idx = subset_idx, getoutvar = FALSE, getXmat = TRUE)
      P1 <- predictP1(self$model.fit,
                      ParentObject = self,
                      DataStorageObject = data,
                      subset_idx = subset_idx,
                      n = self$ParentModel$n
                      # Xmat = private$Xmat,
                      # DataStorageObject = data,
                      )
      return(P1)
    },

    # Sets Xmat, Yvals, evaluates subset and performs correct subseting of data
    # everything is performed using data$ methods (data is of class DataStorageClass)
    setdata = function(data, subset_idx, getoutvar = TRUE, getXmat = TRUE) {
      assert_that(is.DataStorageClass(data))
      if (getoutvar) private$Yvals <- data$get.outvar(subset_idx, self$ParentModel$outvar) # Always a vector
      if (getXmat) self$define.Xmat(data, subset_idx)
      return(invisible(self))
    },

    define.Xmat = function(data, subset_idx) {
      predvars <- self$ParentModel$predvars
      if (sum(subset_idx) == 0L) {  # When nrow(Xmat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
        Xmat <- matrix(, nrow = 0L, ncol = (length(predvars) + 1))
        colnames(Xmat) <- c("Intercept", predvars)
      } else {
        # *** THIS IS THE ONLY LOCATION IN THE PACKAGE WHERE CALL TO DataStorageClass$get.dat.sVar() IS MADE ***
        if (length(predvars)==0L) {
          Xmat <- as.matrix(rep.int(1L, sum(subset_idx)), ncol=1)
        } else {
          Xmat <- as.matrix(cbind(Intercept = 1, data$get.dat.sVar(subset_idx, predvars)))
        }
        colnames(Xmat)[1] <- "Intercept"
        # To find and replace misvals in Xmat:
        if (self$ParentModel$ReplMisVal0) Xmat[gvars$misfun(Xmat)] <- gvars$misXreplace
      }
      private$Xmat <- Xmat
      return(invisible(self))
    }
  ),

  active = list( # 2 types of active bindings (w and wout args)
    emptydata = function() { private$Xmat <- NULL},
    emptyY = function() { private$Yvals <- NULL},
    getXmat = function() {private$Xmat},
    getY = function() {private$Yvals}
  ),

  private = list(
    Xmat = NULL,
    Yvals = NULL
  )
)