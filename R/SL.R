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
