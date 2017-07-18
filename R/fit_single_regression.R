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