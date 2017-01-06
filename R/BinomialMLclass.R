# # Generic for fitting the logistic (binomial family) GLM model
# h2ofit <- function(fit, ...) UseMethod("h2ofit")

# take a list of args, take a function body and return only the args that belong to function signature
keep_only_fun_args <- function(Args, fun) {
  keepArgs <- intersect(names(Args), names(formals(fun))) # captures optional arguments given by user
  if (length(keepArgs) > 0) {
    Args <- Args[keepArgs]
  } else {
    Args <- NULL
  }
  return(Args)
}

# 1. replace any arg in mainArgs if it also appears in userArgs
# 2. add any arg from userArgs that also appears in formals(fun) of function
replace_add_user_args <- function(mainArgs, userArgs, fun) {
  replaceArgs <- intersect(names(mainArgs), names(userArgs)) # captures main arguments that were overridden by user
  if(length(replaceArgs) > 0) {
    mainArgs[replaceArgs] <- userArgs[replaceArgs]
    userArgs[replaceArgs] <- NULL
  }
  newArgs <- intersect(names(formals(fun)), names(userArgs)) # captures any additional args given by user that are not in mainArgs
  if (length(newArgs) > 0) {
    mainArgs <- c(mainArgs, userArgs[newArgs])
  }
  return(mainArgs)
}

# ---------------------------------------------------------------------------
# for running logistic regression with continuous outcome range [0-1]
# ---------------------------------------------------------------------------
# S3 method for fitting h2o GLM with binomial() family (logistic regression):
# use solver="L_BFGS" when doing classification and use "IRLSM" when not
fit.h2oglm <- function(fit.class, fit, training_frame, y, x, model_contrl, ...) {
  h2o::h2o.no_progress()
  mainArgs <- list(x = x, y = y, training_frame = training_frame,
                  intercept = TRUE,
                  family = "binomial",
                  standardize = TRUE,
                  # standardize = FALSE,
                  solver = "L_BFGS",
                  # solver = "IRLSM",
                  # solver = "COORDINATE_DESCENT",
                  # solver = "COORDINATE_DESCENT_NAIVE",
                  lambda = 0L,
                  max_iterations = 100,
                  ignore_const_cols = FALSE,
                  missing_values_handling = "Skip")

  mainArgs <- replace_add_user_args(mainArgs, model_contrl, fun = h2o::h2o.glm)
  model.fit <- do.call(h2o::h2o.glm, mainArgs)
  # assign the fitted coefficients in correct order (same as predictor order in x (former predvars))
  out_coef <- vector(mode = "numeric", length = length(x)+1)
  out_coef[] <- NA
  names(out_coef) <- c("Intercept", x)
  out_coef[names(model.fit@model$coefficients)] <- model.fit@model$coefficients
  fit$coef <- out_coef;
  fit$linkfun <- "logit_linkinv";

  fit$fitfunname <- "h2o.glm";
  confusionMat <- h2o::h2o.confusionMatrix(model.fit)
  fit$nobs <- confusionMat[["0"]][3]+confusionMat[["1"]][3]
  fit$H2O.model.object <- model.fit

  if (gvars$verbose) {
    print("h2oglm fits:")
    print(fit$coef)
  }

  # class(fit) <- c(class(fit)[1], c("glmfit"))
  class(fit) <- c(class(fit)[1], c("H2Omodel"))
  return(fit)
}

# S3 method for h2o randomForest fit (Random Forest):
fit.h2orandomForest <- function(fit.class, fit, training_frame, y, x, model_contrl, ...) {
  h2o::h2o.no_progress()
  mainArgs <- list(x = x, y = y, training_frame = training_frame,
                   ntrees = 100,
                   balance_classes = TRUE,
                   ignore_const_cols = FALSE)

  mainArgs <- replace_add_user_args(mainArgs, model_contrl, fun = h2o::h2o.randomForest)
  model.fit <- do.call(h2o::h2o.randomForest, mainArgs)
  fit$coef <- NULL;
  fit$fitfunname <- "h2o.randomForest";
  confusionMat <- h2o::h2o.confusionMatrix(model.fit)
  fit$nobs <- confusionMat[["0"]][3]+confusionMat[["1"]][3]
  fit$H2O.model.object <- model.fit
  class(fit) <- c(class(fit)[1], c("H2Omodel"))
  return(fit)
}

# S3 method for h2o gbm fit, takes BinDat data object:
# use "bernoulli" when doing classification and use "gaussian" when not
fit.h2ogbm <- function(fit.class, fit, training_frame, y, x, model_contrl, ...) {
  h2o::h2o.no_progress()
  mainArgs <- list(x = x, y = y, training_frame = training_frame,
                   distribution = "bernoulli",
                   # distribution = "gaussian",
                   ntrees = 100,
                   balance_classes = TRUE,
                   ignore_const_cols = FALSE)

  mainArgs <- replace_add_user_args(mainArgs, model_contrl, fun = h2o::h2o.gbm)
  model.fit <- do.call(h2o::h2o.gbm, mainArgs)
  fit$coef <- NULL;
  fit$fitfunname <- "h2o.gbm";
  confusionMat <- h2o::h2o.confusionMatrix(model.fit)
  fit$nobs <- confusionMat[["0"]][3]+confusionMat[["1"]][3]
  fit$H2O.model.object <- model.fit
  class(fit) <- c(class(fit)[1], c("H2Omodel"))
  return(fit)
}

# S3 method for h2o deeplearning fit, takes BinDat data object:
# use "bernoulli" when doing classification and use "gaussian" when doing regression
fit.h2odeeplearning <- function(fit.class, fit, training_frame, y, x, model_contrl, ...) {
  h2o::h2o.no_progress()
  mainArgs <- list(x = x, y = y, training_frame = training_frame,
                   distribution = "bernoulli",
                   # distribution = "gaussian",
                   balance_classes = TRUE,
                   ignore_const_cols = FALSE)

  mainArgs <- replace_add_user_args(mainArgs, model_contrl, fun = h2o::h2o.gbm)
  model.fit <- do.call(h2o::h2o.deeplearning, mainArgs)
  fit$coef <- NULL;
  fit$fitfunname <- "h2o.deeplearning";
  confusionMat <- h2o::h2o.confusionMatrix(model.fit)
  fit$nobs <- confusionMat[["0"]][3]+confusionMat[["1"]][3]
  fit$H2O.model.object <- model.fit
  class(fit) <- c(class(fit)[1], c("H2Omodel"))
  return(fit)
}

# ----------------------------------------------------------------
# Prediction for h2ofit objects, predicts P(A = 1 | newXmat)
# ----------------------------------------------------------------
predictP1.H2Omodel <- function(m.fit, ParentObject, DataStorageObject, subset_idx, n, ...) {
  assert_that(!is.null(subset_idx))
  if (!missing(DataStorageObject)) {
    # rows_subset <- which(subset_idx)
    data <- DataStorageObject
    outvar <- m.fit$params$outvar
    predvars <- m.fit$params$predvars

    # 1. works on a single core, but fails in parallel:
    subsetH2Oframe <- data$fast.load.to.H2O(data$dat.sVar[subset_idx, c(outvar, predvars), with = FALSE],
                                            saveH2O = FALSE,
                                            destination_frame = "subsetH2Oframe")
    #2. old, slower approach, but may work on many cores (since data is loaded only once)
    # subsetH2Oframe <- data$H2O.dat.sVar[subset_idx, c(outvar, predvars)]
    # old version of setting data, no longer used:
    # ParentObject$setdata(data, subset_idx = subset_idx, getoutvar = FALSE, getXmat = FALSE)

  } else {
    subsetH2Oframe <- ParentObject$getsubsetH2Oframe
  }

  pAout <- rep.int(gvars$misval, n)
  if (length(subset_idx) > 0) {
  # if (sum(subset_idx) > 0) {
    predictFrame <- h2o::h2o.predict(m.fit$H2O.model.object, newdata = subsetH2Oframe)
    # predictObject <- predict(m.fit$H2O.model.object, newdata = subsetH2Oframe)
    # predictFrame <- predictObject$pred
    if ("p1" %in% colnames(predictFrame)) {
      pAout[subset_idx] <- as.vector(predictFrame[,"p1"])
    } else {
      pAout[subset_idx] <- as.vector(predictFrame[,"predict"])
    }
  }
  return(pAout)
}

# ----------------------------------------------------------------
# Prediction for h2ofit objects, predicts P(A = 1 | newXmat)
# ----------------------------------------------------------------
predictP1.H2Oensemblemodel <- function(m.fit, ParentObject, DataStorageObject, subset_idx, n, ...) {
  assert_that(!is.null(subset_idx))
  if (!missing(DataStorageObject)) {
    # rows_subset <- which(subset_idx)
    data <- DataStorageObject
    outvar <- m.fit$params$outvar
    predvars <- m.fit$params$predvars

    # 1. works on a single core, but fails in parallel:
    subsetH2Oframe <- data$fast.load.to.H2O(data$dat.sVar[subset_idx, c(outvar, predvars), with = FALSE],
                                            saveH2O = FALSE,
                                            destination_frame = "subsetH2Oframe")
    #2. old, slower approach, but may work on many cores (since data is loaded only once)
    # subsetH2Oframe <- data$H2O.dat.sVar[subset_idx, c(outvar, predvars)]
    # old version of setting data, no longer used:
    # ParentObject$setdata(data, subset_idx = subset_idx, getoutvar = FALSE, getXmat = FALSE)

  } else {
    subsetH2Oframe <- ParentObject$getsubsetH2Oframe
  }

  pAout <- rep.int(gvars$misval, n)
  # if (sum(subset_idx) > 0) {
  if (length(subset_idx) > 0) {
    predictObject <- predict(m.fit$H2O.model.object, newdata = subsetH2Oframe)
    predictFrame <- predictObject$pred
    if ("p1" %in% colnames(predictFrame)) {
      pAout[subset_idx] <- as.vector(predictFrame[,"p1"])
    } else {
      pAout[subset_idx] <- as.vector(predictFrame[,"predict"])
    }
  }
  return(pAout)
}

# IMPLEMENTING NEW CLASS FOR BINARY REGRESSION THAT USES h2o
# NEEDS TO be able to pass on THE REGRESSION SETTINGS FOR h2o-specific functions
BinomialH2O  <- R6Class(classname = "BinomialH2O",
  inherit = BinomialGLM,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    fit.class = character(),
    model.fit = list(coef = NA, fitfunname = NA, linkfun = NA, nobs = NA, params = NA, H2O.model.object = NA),
    outfactors = NA,
    nfolds = 5,

    initialize = function(fit.algorithm, fit.package, ParentModel, ...) {
      self$ParentModel <- ParentModel
      self$classify <- ParentModel$classify
      self$model_contrl <- ParentModel$model_contrl
      assert_that("h2o" %in% fit.package)

      # val <- checkpkgs(pkgs = c("h2o"))
      if (!requireNamespace("h2o", quietly = TRUE)) {
        stop(
"Package h2o is needed for machine learning. Please make sure to check that Java is also installed.
To install the CRAN version of h2o package (might not be the latest), type this into R terminal:
  install.packages(\"h2o\")
................................
For instructions on using h2o and installing its latest version for R go to:
  https://s3.amazonaws.com/h2o-release/h2o/master/latest.html (click on 'INSTALL IN R')",
        call. = FALSE)
      }

      self$fit.class <- fit.algorithm
      class(self$fit.class) <- c(class(self$fit.class), "h2o" %+% self$fit.class)
      invisible(self)
    },

    fit = function(data, outvar, predvars, subset_idx, ...) {
      assert_that(is.DataStorageClass(data))
      self$setdata(data, subset_idx, self$classify, ...)
      if ((length(predvars) == 0L) || (length(subset_idx) == 0L) || (length(self$outfactors) < 2L)) {
      # if ((length(predvars) == 0L) || (sum(subset_idx) == 0L) || (length(self$outfactors) < 2L)) {
        message("unable to run " %+% self$fit.class %+% " with h2o for: intercept only models or designmat with zero rows or  constant outcome (y) ...")
        class(self$model.fit) <- "try-error"
        self$emptydata
        self$emptyY
        return(self$model.fit)
      }

      self$model.fit$params <- self$params
      self$model.fit <- try(
                    fit(self$fit.class, self$model.fit,
                           training_frame = private$subsetH2Oframe,
                           y = outvar, x = predvars,
                           model_contrl = self$model_contrl, fold_column = data$fold_column, ...),
                    silent = FALSE)

      if (inherits(self$model.fit, "try-error")) {
        self$emptydata
        self$emptyY
        return(self$model.fit)
      }
      return(self$model.fit)
    },

    setdata = function(data, subset_idx, classify = TRUE, getoutvar = TRUE, ...) {
      outvar <- self$ParentModel$outvar
      predvars <- self$ParentModel$predvars
      # rows_subset <- which(subset_idx)
      # a penalty for being able to obtain predictions from predictAeqA() right after fitting: need to store Yvals
      if (getoutvar) private$Yvals <- data$get.outvar(subset_idx, outvar) # Always a vector

      if (self$fit.class %in% "SuperLearner") {
        if (!is.null(self$model_contrl$nfolds)) {
          self$nfolds <- as.integer(self$model_contrl$nfolds)
        } else {
          self$model_contrl$nfolds <- self$nfolds
        }
        data$define_CVfolds(nfolds = self$nfolds, fold_column = "fold_id", seed = self$model_contrl$seed)
      }
      # ---------------------------------------------------
      # ***!!!!NEED TO ALLOW USING BOTH VERSIONS BELOW!!!!***
      # ---------------------------------------------------
      # 1. works on single core but fails in parallel:
      load_subset_t <- system.time(
        subsetH2Oframe <- data$fast.load.to.H2O(data$dat.sVar[subset_idx, c(outvar, predvars, data$fold_column), with = FALSE],
                                                saveH2O = FALSE,
                                                destination_frame = "newH2Osubset")
      )
      if (gvars$verbose) {
        print("time to subset and load data into H2OFRAME: "); print(load_subset_t)
      }
      # 2. old version, less efficient when subsetting really large frames:
      # ---------------------------------------------------
      # everytime you do a subset H2O creates a frame on its cluster
      # if you run this in parallel, it will try to assign the same name to two different subsets calls below
      # resuling in an error
      # ---------------------------------------------------
      # subset_t <- system.time(
      #   subsetH2Oframe <- data$H2O.dat.sVar[subset_idx, c(outvar, predvars)]
      # )
      # if (gvars$verbose) {
      #   print("time to subset data into H2OFRAME: "); print(subset_t)
      # }

      # **** WILL CREATE A TEMPORARY h2o FRAME ****
      self$outfactors <- as.vector(h2o::h2o.unique(subsetH2Oframe[, outvar]))
      # Below being TRUE implies that the conversion to H2O.FRAME produced errors, since there should be no NAs in the source subset data
      if (any(is.na(self$outfactors))) stop("FOUND NA OUTCOMES IN H2OFRAME WHEN THERE WERE NOT SUPPOSED TO BE ANY")

      if (classify) {
        if (length(self$outfactors) > 2L) stop("cannot run binary regression/classification for outcome with more than 2 categories")
        subsetH2Oframe[, outvar] <- h2o::as.factor(subsetH2Oframe[, outvar])
      }
      private$subsetH2Oframe <- subsetH2Oframe
      return(invisible(self))
    }
  ),

  active = list( # 2 types of active bindings (w and wout args)
    emptydata = function() { private$subsetH2Oframe <- NULL},
    getsubsetH2Oframe = function() {private$subsetH2Oframe}
  ),

  private = list(
    subsetH2Oframe = NULL
  )
)
