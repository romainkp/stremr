# Generic for fitting the logistic (binomial family) GLM model
h2ofit <- function(fit, ...) UseMethod("h2ofit")

h2ofit.h2oSL <- function(fit, DataStorageObject, outvar, predvars, subset_idx, ...) {
  # ...
}

h2ofit.h2oGLM <- function(fit, DataStorageObject, outvar, predvars, subset_idx, ...) {
  h2o.glmfit <- glmfit.h2oglm(fit, DataStorageObject, outvar, predvars, subset_idx, ...)
  return(h2o.glmfit)
}

# S3 method for h2o RFs fit, takes BinDat data object:
h2ofit.h2oRF <- function(fit, DataStorageObject, outvar, predvars, subset_idx, ...) {
  if (gvars$verbose) print("calling h2o.randomForest...")
  # yname <- BinDatObject$outvar
  # xnames <- BinDatObject$predvars
  # subset_idx <- which(BinDatObject$subset_idx)
  # if (length(subset_idx) == 0L) { # Xmat has 0 rows: return NA`s and avoid throwing exception
  #   model.fit <- list(coef = rep.int(NA_real_, length(xnames)))
  # } else if (length(xnames) == 0L) {
  #   return(logisfit.speedglm(BinDatObject))
  # } else {
  # Random Forests:
  model.fit <- h2o::h2o.randomForest(x = predvars,
                                     y = outvar,
                                     training_frame = DataStorageObject$H2O.dat.sVar[subset_idx,],
                                     ntree = 100)

  fit$coef <- model.fit$coef; fit$fitfunname <- "h2o.randomForest"; fit$nobs <- length(subset_idx);
  fit$H2O.model.object <- model.fit
  # fit <- list(coef = NA, linkfun = NA, fitfunname = "h2o.randomForest", nobs = length(subset_idx), H2O.model.object = model.fit)
  # if (gvars$verbose) print(fit$coef)
  class(fit) <- c(class(fit)[1], c("h2ofit"))
  return(fit)
}

# S3 method for h2o GBM fit, takes BinDat data object:
h2ofit.h2oGBM <- function(fit, DataStorageObject, outvar, predvars, subset_idx, ...) {
  if (gvars$verbose) print("calling h2o.gbm...")
    # yname <- BinDatObject$outvar
    # xnames <- BinDatObject$predvars
    # subset_idx <- which(BinDatObject$subset_idx)
    # if (length(subset_idx) == 0L) { # Xmat has 0 rows: return NA`s and avoid throwing exception
    #   model.fit <- list(coef = rep.int(NA_real_, length(xnames)))
    # } else if (length(xnames) == 0L) {
    #   return(logisfit.speedglmS3(BinDatObject))
    # } else {
    # for GBM to run need to make outcome into a factor:
    newH2Oframe <- DataStorageObject$H2O.dat.sVar[subset_idx,]
    newH2Oframe[,yname] <- h2o::as.factor(newH2Oframe[,yname])
    model.fit <- h2o::h2o.gbm(x = predvars,
                          y = outvar,
                          training_frame = newH2Oframe,
                          distribution = "bernoulli")

    fit$coef <- model.fit$coef; fit$fitfunname <- "h2o.gbm"; fit$nobs <- length(subset_idx); fit$H2O.model.object <- model.fit
    # fit <- list(coef = NA, linkfun = NA, fitfunname = "h2o.gbm", nobs = length(subset_idx), H2O.model.object = model.fit)
    # if (gvars$verbose) print(fit$coef)
    class(fit) <- c(class(fit)[1], c("h2ofit"))
    return(fit)
}


# IMPLEMENTING NEW CLASS FOR BINARY REGRESSION THAT DIRECTLY USES PRELOADED H20FRAME
# DOES NOT NEED TO RETRIEVE ANY COVARS, NEEDS TO ONLY EVALUATE THE SUBSET_IDX
# NEEDS TO be able to pass on THE REGRESSION SETTINGS FOR h2o-specific functions
# POSSIBLY NEEDS A SEPARATE RegressionClass
BinomialH2O  <- R6Class(classname = "BinomialH2O",
  inherit = BinomialGLM,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(

    # TO DO: THIS WILL CONTAIN ADDITIONAL USER-SPEC'ED CONTROLS/ARGS PASSED ON TO h2o or h2oEnsemble
    model.controls = NULL,
    model.fit = list(coef = NA, fitfunname = NA, nobs = NA, H2O.model.object = NA),

    fit.class = c("h2oSL", "h2oGLM", "h2oRF", "h2oGBM"),

    initialize = function(fit.algorithm, fit.package, ParentModel, ...) {
      self$ParentModel <- ParentModel
      assert_that(any(c("h2o", "h2oglm") %in% fit.package))
      self$fit.class <- fit.algorithm
      class(self$model.fit) <- c(class(self$model.fit), fit.algorithm)

      invisible(self)
    },

    fit = function(data, outvar, predvars, subset_idx, ...) {
      self$setdata(data, subset_idx = subset_idx, getoutvar = FALSE, getXmat = FALSE)
      # Xmat has 0 rows: return NA's and avoid throwing exception:
      model.fit <- try(
                      h2ofit(self$model.fit,
                               DataStorageObject = data,
                               outvar = outvar,
                               predvars = predvars,
                               subset_idx = subset_idx, ...),
                      silent = TRUE)
      if (inherits(model.fit, "try-error")) { # failed, need to define the Xmat now and try fitting speedglm/glm
        class(self$model.fit)[2] <- "speedglm"
        warning(self$fit.class %+% " failed, falling back on speedglm...")
        # self$setdata(data, subset_idx = subset_idx, getoutvar = FALSE, getXmat = TRUE)
        model.fit <- super$fit(data, outvar, predvars, subset_idx, ...)
      }

      self$model.fit <- model.fit
      return(self$model.fit)
    },

    predictP1 = function(data, subset_idx) {
      # browser()
      # self$setdata(data, subset_idx = subset_idx, getoutvar = FALSE, getXmat = FALSE)
      P1 <- predictP1(self$model.fit,
                      ParentObject = self,
                      DataStorageObject = data,
                      subset_idx = subset_idx,
                      n = self$ParentModel$n)
      return(P1)
    }

    # # Sets Xmat, Yvals, evaluates subset and performs correct subseting of data
    # # everything is performed using data$ methods (data is of class DataStorageClass)
    # setdata = function(data, subset_idx, getoutvar = TRUE,  ...) {
    #   assert_that(is.DataStorageClass(data))
    #   if (getoutvar) private$Yvals <- data$get.outvar(subset_idx, self$ParentModel$outvar) # Always a vector
    #   # NO NEED FOR Xmat with H2O, but might need to do some manipulations on data (load biniried matrix in H2O.FRAME)
    #   # self$define.Xmat(data, subset_idx)
    #   return(invisible(self))
    # }
  )
)







