# Generic for fitting the logistic (binomial family) GLM model
h2ofit <- function(fit, ...) UseMethod("h2ofit")

# S3 method for fitting h2o GLM with binomial() family (logistic regression):
h2ofit.h2oGLM <- function(fit, DataStorageObject, subsetH2Oframe, outvar, predvars, rows_subset, ...) {
  if (gvars$verbose) print("calling h2o.glm...")
  model.fit <- h2o::h2o.glm(x = predvars,
                            y = outvar,
                            intercept = TRUE,
                            training_frame = subsetH2Oframe,
                            family = "binomial",
                            standardize = TRUE,
                            solver = c("L_BFGS"),
                            # solver = c("IRLSM"),
                            # remove_collinear_columns = TRUE,
                            max_iterations = 50,
                            lambda = 0L)

    # print("h2o.glm fit"); print(model.fit)
    # print(model.fit@parameters)
    # print(model.fit@allparameters)
    # str(model.fit@model)

  # assign the fitted coefficients in correct order (same as predictor order in predvars)
  out_coef <- vector(mode = "numeric", length = length(predvars)+1)
  out_coef[] <- NA
  names(out_coef) <- c("Intercept", predvars)
  out_coef[names(model.fit@model$coefficients)] <- model.fit@model$coefficients

  fit$coef <- out_coef;
  fit$fitfunname <- "h2o.glm";
  fit$linkfun <- "logit_linkinv";
  fit$nobs <- length(rows_subset);
  fit$H2O.model.object <- model.fit

  if (gvars$verbose) {
    print("h2oglm fits:")
    print(fit$coef)
  }

  class(fit) <- c(class(fit)[1], c("glmfit"))
  # class(fit) <- c(class(fit)[1], c("h2ofit"))
  return(fit)
}

# S3 method for h2o RFs fit (Random Forest):
h2ofit.h2oRF <- function(fit, DataStorageObject, subsetH2Oframe, outvar, predvars, rows_subset, ...) {
  if (gvars$verbose) print("calling h2o.randomForest...")
  model.fit <- h2o::h2o.randomForest(x = predvars,
                                     y = outvar,
                                     training_frame = subsetH2Oframe,
                                     ntree = 100,
                                     ignore_const_cols = FALSE)
  fit$coef <- NA;
  fit$fitfunname <- "h2o.randomForest";
  fit$nobs <- length(rows_subset);
  fit$H2O.model.object <- model.fit

  class(fit) <- c(class(fit)[1], c("h2ofit"))
  return(fit)
}

# S3 method for h2o GBM fit, takes BinDat data object:
h2ofit.h2oGBM <- function(fit, DataStorageObject, subsetH2Oframe, outvar, predvars, rows_subset, ...) {
  if (gvars$verbose) print("calling h2o.gbm...")
    model.fit <- h2o::h2o.gbm(x = predvars,
                              y = outvar,
                              training_frame = subsetH2Oframe,
                              distribution = "bernoulli",
                              ignore_const_cols = FALSE)
    fit$coef <- NA;
    fit$fitfunname <- "h2o.gbm";
    fit$nobs <- length(rows_subset);
    fit$H2O.model.object <- model.fit
    class(fit) <- c(class(fit)[1], c("h2ofit"))
    return(fit)
}

h2ofit.h2oSL <- function(fit, DataStorageObject, outvar, predvars, subset_idx, ...) {
  # ...
  # ... SuperLearner TO BE IMPLEMENTED ...
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
      model.fit <- self$model.fit

      if ((length(predvars) == 0L) || (sum(subset_idx) == 0L)) {
        class(model.fit) <- "try-error"
      } else {
        rows_subset <- which(subset_idx)
        subsetH2Oframe <- data$H2O.dat.sVar[rows_subset, ]

        outfactors <- as.vector(h2o::h2o.unique(subsetH2Oframe[, outvar]))
        # Below being TRUE implies that the conversion to H2O.FRAME produced errors, since there are no NAs in the original data
        NAfactors <- any(is.na(outfactors))

        # fixing bug in h2o frame
        if (NAfactors) {
          NA_idx_h2o <- which(as.logical(is.na(subsetH2Oframe[,outvar])))
          orig.vals <- data$dat.sVar[rows_subset, ][NA_idx_h2o, outvar, with = FALSE][[outvar]]
          # require("h2o")
          subsetH2Oframe[NA_idx_h2o, outvar] <- orig.vals[1]
          outfactors <- as.vector(h2o::h2o.unique(subsetH2Oframe[, outvar]))
          NAfactors <- any(is.na(outfactors))
        }

        if (length(outfactors) < 2L | NAfactors) {
          class(model.fit) <- "try-error"
        } else if (length(outfactors) > 2L) {
          stop("Attempting to run Binary regression/classification for outcome with more than 2 categories")
        }
      }

      if (!inherits(model.fit, "try-error")) {
        subsetH2Oframe[, outvar] <- h2o::as.factor(subsetH2Oframe[, outvar])
        model.fit <- try(
                h2ofit(self$model.fit,
                       DataStorageObject = data,
                       subsetH2Oframe = subsetH2Oframe,
                       outvar = outvar,
                       predvars = predvars,
                       rows_subset = rows_subset, ...),
                silent = TRUE)
      }

      if (inherits(model.fit, "try-error")) { # failed, need to define the Xmat now and try fitting speedglm/glm
        message(self$fit.class %+% " failed, falling back on speedglm...")
        class(self$model.fit)[2] <- "speedglm"
        model.fit <- super$fit(data, outvar, predvars, subset_idx, ...)
      }

      self$model.fit <- model.fit
      return(self$model.fit)
    }

  )
)







