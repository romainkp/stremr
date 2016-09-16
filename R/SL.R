#' h2o glm wrapper for non-negative least squares
#'
#' This wrapper is used as a default metalearner for SuperLearner / h2oEnsemble.
#' @param ... All arguments which will be passed down to \code{h2oEnsemble::h2o.glm.wrapper(...)}
#' @param non_negative Flag enabling non-negative least squares
#' @return ...
#' @export
h2o.glm_nn <- function(..., non_negative = TRUE) h2oEnsemble::h2o.glm.wrapper(..., non_negative = non_negative)

# Train a model using a single h2o learner (user-spec) with cross-validation & keep CV predictions
SLfit.h2oLearner <- function(learner, training_frame, y, x, family = "binomial", fold_column, model_contrl, ...) {
  # if (is.numeric(seed)) set.seed(seed)  #If seed given, set seed prior to next step
  h2o::h2o.no_progress()

  if (("x" %in% names(formals(learner))) && (as.character(formals(learner)$x)[1] != "")) {
    # Special case where we pass a subset of the colnames, x, in a custom learner function wrapper
    model.fit <- match.fun(learner)(y = y, training_frame = training_frame, validation_frame = NULL, family = family, fold_column = fold_column, keep_cross_validation_folds = TRUE)
  } else {
    # Use all predictors in training set for training
    model.fit <- match.fun(learner)(y = y, x = x, training_frame = training_frame, validation_frame = NULL, family = family, fold_column = fold_column, keep_cross_validation_folds = TRUE)
  }

  fit <- vector(mode = "list")
  fit$fitfunname <- learner;
  if (gvars$verbose) {
    print("grid search fitted models:"); print(model.fit)
  }
  fit$H2O.model.object <- model.fit
  class(fit) <- c(class(fit)[1], c("H2Omodel"))
  return(fit)
}

# S3 method for h2o deeplearning fit, takes BinDat data object:
# use "bernoulli" when doing classification and use "gaussian" when doing regression
SLfit.h2ogrid <- function(grid.algorithm, training_frame, y, x, family = "binomial", fold_column, model_contrl, ...) {
  h2o::h2o.no_progress()
  mainArgs <- list(x = x, y = y, training_frame = training_frame,
                  intercept = TRUE,
                  seed = 1,
                  fold_column = fold_column,
                  # fold_assignment = "Modulo",
                  keep_cross_validation_predictions = TRUE,
                  family = family,
                  standardize = TRUE,
                  solver = "L_BFGS",
                  lambda = 0L,
                  max_iterations = 100,
                  ignore_const_cols = FALSE,
                  missing_values_handling = "Skip")

  if (is.null(grid.algorithm)) stop("must specify 'grid.algorithm' name when running 'h2o.grid'")
  if (!is.character(grid.algorithm)) stop("'grid.algorithm' must be a string naming the grid.algorithm for 'h2o.grid'")
  algo_fun_name <- "h2o."%+%grid.algorithm
  if (!exists(algo_fun_name)) stop("could not locate the function %+% " %+% grid.algorithm)

  algo_fun <- get0(algo_fun_name, mode = "function", inherits = TRUE)
  mainArgs <- keep_only_fun_args(mainArgs, fun = algo_fun)   # Keep only the relevant args in mainArgs list:

  mainArgs <- replace_add_user_args(mainArgs, model_contrl, fun = algo_fun) # Add user args that pertain to this specific learner:
  mainArgs$algorithm <- grid.algorithm
  mainArgs$search_criteria <- model_contrl[["search_criteria"]]
  mainArgs$hyper_params <- model_contrl[[grid.algorithm]]

  if (is.null(mainArgs$hyper_params)) stop("must specify hyper parameters for grid search with '" %+% algo_fun_name %+% "' by defining a SuperLearner params list item named '" %+% grid.algorithm %+% "'")

  if (!is.null(mainArgs$hyper_params[["search_criteria"]])) {
    mainArgs$search_criteria <- mainArgs$hyper_params[["search_criteria"]]
    mainArgs$hyper_params[["search_criteria"]] <- NULL
  }
  if (is.null(mainArgs$search_criteria)) stop("must specify 'search_criteria' when running 'h2o.grid' for grid.algorithm " %+% grid.algorithm)

  # Remove any args from mainArgs that also appear in hyper_params:
  common_hyper_args <- intersect(names(mainArgs), names(mainArgs$hyper_params))
  if(length(common_hyper_args) > 0) mainArgs <- mainArgs[!(names(mainArgs) %in% common_hyper_args)]

  if (gvars$verbose) {
    print("running h2o.grid grid.algorithm: "); print(grid.algorithm)
  }

  model.fit <- do.call(h2o::h2o.grid, mainArgs)

  fit <- vector(mode = "list")
  fit$fitfunname <- "h2o.h2ogrid";
  if (gvars$verbose) {
    print("grid search fitted models:"); print(model.fit)
  }
  fit$H2O.model.object <- model.fit
  class(fit) <- c(class(fit)[1], c("H2Omodel"))
  return(fit)
}

fit.h2oSuperLearner <- function(fit.class, fit, training_frame, y, x, model_contrl, fold_column, ...) {
  val <- checkpkgs(pkgs = c("h2oEnsemble"))
  if (is.null(fold_column)) stop("must define the column of CV fold IDs using data$define_CVfolds()")

  if (!is.null(model_contrl$load.ensemble) && model_contrl$load.ensemble) {
    if (is.null(model_contrl$ensemble.dir.path) || (model_contrl$ensemble.dir.path %in% "")) {
      stop("when loading ensemble must specify the directory path with 'ensemble.dir.path' parameter")
    }
    stacked.fit <- h2oEnsemble::h2o.load_ensemble(path = model_contrl$ensemble.dir.path, import_levelone = FALSE)
    # stacked.fit <- h2oEnsemble::h2o.load_ensemble(path = model_contrl$ensemble.dir.path, import_levelone = TRUE)

  } else {
    family <- model_contrl$family
    grid.algorithms <- model_contrl$grid.algorithm
    learners <- model_contrl$learner
    if (is.null(family)) family <- "binomial"
    # Will put all fitted models in a single list for stacking
    fitted_models_all <- NULL
    nfolds <- model_contrl$nfolds
    model_contrl$nfolds <- NULL

    if (is.null(grid.algorithms) && is.null(learners)) {
      stop("must specify either 'grid.algorithm' or 'learner' when performing estimation with SuperLearner")
    }

    if (!is.null(grid.algorithms)) {
      if (!is.character(grid.algorithms)) stop("'grid.algorithm' must be a vector of strings naming the grid.algorithms to use in 'h2o.grid'")
      fitted_models <- vector(mode = "list", length = length(grid.algorithms))
      names(fitted_models) <- grid.algorithms
      for (grid.algorithm in grid.algorithms) {
        grid_model_fit <- SLfit.h2ogrid(grid.algorithm = grid.algorithm, training_frame = training_frame, y = y, x = x, family = family, fold_column = fold_column, model_contrl = model_contrl, ...)
        grid_model_H2O <- grid_model_fit$H2O.model.object
        fitted_models[[grid.algorithm]] <- lapply(grid_model_H2O@model_ids, function(model_id) h2o::h2o.getModel(model_id))
      }
      for (grid.algorithm in grid.algorithms) {
        fitted_models_all <- c(fitted_models_all, fitted_models[[grid.algorithm]])
      }
    }

    if (!is.null(learners)) {
      if (!is.character(learners)) stop("'learner' must be a vector of strings naming specific wrappers/learners")
      fitted_models_l <- vector(mode = "list", length = length(learners))
      names(fitted_models_l) <- learners
      for (learner in learners) {
        learner_fit <- SLfit.h2oLearner(learner = learner, training_frame = training_frame, y = y, x = x, family = family, fold_column = fold_column, model_contrl = model_contrl, ...)
        learner_model_H2O <- learner_fit$H2O.model.object
        fitted_models_l[[learner]] <- learner_model_H2O
      }
      fitted_models_all <- c(fitted_models_all, fitted_models_l)
    }

    # to by-pass error check in h2o.stack:
    for (idx in seq_along(fitted_models_all)) {
      fitted_models_all[[idx]]@allparameters$fold_assignment <- "Modulo"
      fitted_models_all[[idx]]@allparameters$nfolds <- nfolds
    }

    # Specify a defalt GLM as the metalearner
    metalearner <- model_contrl$metalearner
    if (is.null(metalearner)) metalearner <- "h2o.glm_nn"
    stacked.fit <- h2oEnsemble::h2o.stack(models = fitted_models_all, response_frame = training_frame[,y], metalearner = metalearner)

    # ----------------------------------------------------------------------------------------------------
    # Saving the fits:
    # ----------------------------------------------------------------------------------------------------
    if (!is.null(model_contrl$save.ensemble) && model_contrl$save.ensemble) {
      if (is.null(model_contrl$ensemble.dir.path) || (model_contrl$ensemble.dir.path %in% "")) {
        stop("when saving ensemble must specify the directory path with 'ensemble.dir.path' parameter")
      }
      h2oEnsemble::h2o.save_ensemble(stacked.fit, path = model_contrl$ensemble.dir.path, force = TRUE)
      # h2oEnsemble::h2o.save_ensemble(stacked.fit, path = model_contrl$ensemble.dir.path, force = TRUE, export_levelone = TRUE)
    }
  }

  # ----------------------------------------------------------------------------------------------------
  # Compute the final SL performance on the training set:
  # ----------------------------------------------------------------------------------------------------
  print("SuperLearner fit:"); print(stacked.fit$metafit)
  perf <- h2oEnsemble::h2o.ensemble_performance(stacked.fit, newdata = training_frame, score_base_models = FALSE)
  print("SuperLearner overall performance (AUC) on the training set: "); print(perf)
  print("SuperLearner overall performance (MSE) on the training set: "); print(perf$ensemble@metrics$MSE)
  # h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)
  # stacked.fit3 <- h2o.metalearn(stacked.fit, metalearner = "h2o.glm_nn")
  # perf3 <- h2o.ensemble_performance(stacked.fit3, newdata = training_frame, score_base_models = FALSE)
  # print(perf3)

  # out_coef <- vector(mode = "numeric", length = length(fitted_models_all)+1)
  out_coef <- vector(mode = "numeric", length = length(stacked.fit$learner)+1)
  out_coef[] <- NA
  names(out_coef) <- names(stacked.fit$metafit@model$coefficients)
  out_coef[names(stacked.fit$metafit@model$coefficients)] <- stacked.fit$metafit@model$coefficients
  names(out_coef)[which(!names(stacked.fit$learner) %in% "")+1] <- names(stacked.fit$learner)[!names(stacked.fit$learner) %in% ""]

  fit$coef <- out_coef;
  fit$linkfun <- NA
  fit$nobs <- nrow(training_frame)

  if (gvars$verbose) {
    print("SuperLearner fits:")
    print(fit$coef)
  }

  fit$fitfunname <- "h2oEnsemble::h2o.stack";
  fit$H2O.model.object <- stacked.fit

  # TO DIRECTLY SAVE ALL MODEL FITS FROM GRID SEARCH (base-learners)
  # fit$fitted_models_all <- fitted_models_all

  class(fit) <- c(class(fit)[1], c("H2Oensemblemodel"))
  return(fit)
}