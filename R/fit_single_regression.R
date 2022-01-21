## ----------------------------------------------------------------------------------
## Call \code{gridisl} and fit a single regression model.
## All model fitting in stremr is performed via this function.
## ----------------------------------------------------------------------------------
fit_single_regression <- function(data, nodes, models, model_contrl, predvars, outvar, subset_idx, outcome_type = "continuous", ...) {

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

  if (is(models, "Lrnr_base")) {
    ## todo: implement an efficient converter from DataStorageClass into sl3_task object (avoiding a copy unless absolutely necessary)
    ## todo: need to somehow transfer the fold column to task, but only if it actually exists
    if (!is.null(data$fold_column)) {
      folds <- data$make_origami_fold_from_column(subset_idx)
    } else {
      folds <- NULL
    }

    task <- sl3::sl3_Task$new(data$dat.sVar[subset_idx, ],
                              covariates = predvars,
                              outcome = outvar,
                              id = data$nodes$IDnode,
                              folds = folds,
                              outcome_type = outcome_type
                              )
    model.fit <- try({models$train(task)})

    if (inherits(model.fit, "try-error") || inherits(model.fit$fit_object, "try-error")) {
      cat("\nsl3 error debugging info; NA is assigned to all predictions:\n");
      print(model.fit)
    }
    # try({
    #   internal_ref <- model.fit$training_task$data
    #   data.table::set(internal_ref, j=names(internal_ref), value=NULL)
    # }, silent = TRUE)

  } else {
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
  }

  if (inherits(model.fit, "try-error") || inherits(model.fit$fit_object, "try-error")) {
    message("...trying to run Lrnr_glm_fast as a backup...")
    task <- sl3::sl3_Task$new(data$dat.sVar[subset_idx, ], covariates = predvars, outcome = outvar, outcome_type = "continuous")
    lrn_model <- sl3::Lrnr_glm$new()
    model.fit <- try(lrn_model$train(task))
    if (inherits(model.fit, "try-error")) {
      cat("\nsl3 error debugging info:\n");
      print(model.fit)
    }
    # try({
    #   internal_ref <- model.fit$training_task$data
    #   data.table::set(internal_ref, j=names(internal_ref), value=NULL)
    # }, silent = TRUE)
  }

  return(model.fit)
}