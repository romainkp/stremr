#----------------------------------------------------------------------------------
# S3 classes for printing model summaries
#----------------------------------------------------------------------------------
#' @importFrom pander pander
NULL

#' Pander method for H2OBinomialMetrics S4 class
#'
#' Prints a H2OBinomialMetrics object in Pandoc's markdown.
#' @param H2OBinomialMetricsObject H2OBinomialMetrics object
#' @param type Character name specifying the type of metric (e.g., "training", "validation", "cross-validation")
#' @return By default this function outputs (see: \code{?cat}) the result.
#' If you would want to catch the result instead, then call \code{pander_return} instead.
#' @export
pander.H2OBinomialMetrics <- function(H2OBinomialMetricsObject, type) {

  h2o_metrics <- H2OBinomialMetricsObject@metrics
  modelID <- h2o_metrics$model$name
  categor <- h2o_metrics$model_category
  # metricsDF <- t(data.frame(h2o_metrics[names(h2o_metrics) %in% c("nobs", "MSE", "RMSE", "mean_residual_deviance", "mae", "rmsle", "logloss", "AUC", "Gini")]))
  metricsDF <- t(data.table::as.data.table(
    h2o_metrics[
      names(h2o_metrics) %in% c("nobs", "MSE", "RMSE", "r2", "logloss", "AUC", "mean_residual_deviance", "mae", "rmsle", "Gini")]
      ))

  metricsDF <- data.table::data.table(
    metric = rownames(metricsDF),
    value = as.numeric(metricsDF[,1]))

  if (H2OBinomialMetricsObject@algorithm == "glm") {
    glmDF <- data.table::data.table(
      metric =
        c("Null Deviance:      ",
          "Residual Deviance:  ",
          "AIC:                "),
      value = c(h2o_metrics$null_deviance,
            h2o_metrics$residual_deviance,
            h2o_metrics$AIC))
    metricsDF <- rbind(metricsDF, glmDF)
  }
  # colnames(metricsDF) <- NULL
  pander::pander(metricsDF, justify = c('left', 'center'), caption = type %+% " data metrics" %+% "; Category: " %+% categor)

  max_matrics <- h2o_metrics$max_criteria_and_metric_scores
  caption <- attributes(max_matrics)$header %+% ": " %+% attributes(max_matrics)$description
  pander::pander(max_matrics, caption = caption %+% " (Model ID: " %+% modelID %+%")")
  # attributes(max_matrics)$formats
  # cat("\nGains/Lift Table: Extract with `h2o.gainsLift(<model>, <data>)` or `h2o.gainsLift(<model>, valid=<T/F>, xval=<T/F>)`")
  if (!is.null(h2o_metrics$gains_lift_table)) {
   # print(h2o_metrics$gains_lift_table)
   gain_tab <- h2o_metrics$gains_lift_table
   pander::pander(gain_tab, caption = attributes(gain_tab)$header %+% " (Model ID: " %+% modelID %+%")")
  }
  return(invisible(H2OBinomialMetricsObject))
}

#' Pander method for H2ORegressionMetrics S4 class
#'
#' Prints a H2ORegressionMetrics object in Pandoc's markdown.
#' @param H2ORegressionMetricsObject H2ORegressionMetrics object
#' @param type Character name specifying the type of metric (e.g., "training", "validation", "cross-validation")
#' @return By default this function outputs (see: \code{?cat}) the result.
#' If you would want to catch the result instead, then call \code{pander_return} instead.
#' @export
pander.H2ORegressionMetrics <- function(H2ORegressionMetricsObject, type) {
  pander.H2OBinomialMetrics(H2ORegressionMetricsObject, type)
}

#' Pander method for H2OGrid S4 class
#'
#' Prints a H2OGrid object in Pandoc's markdown.
#' @param H2OGridObject H2OGrid S4 object
#' @return By default this function outputs (see: \code{?cat}) the result.
#' If you would want to catch the result instead, then call \code{pander_return} instead.
#' @export
pander.H2OGrid <- function(H2OGridObject) {
  grid_summary_tab <- H2OGridObject@summary_table
  caption <- "H2O Grid. " %+% attr(grid_summary_tab, "header") %+% ": " %+% attr(grid_summary_tab, "description")
  pander::pander(data.table::as.data.table(grid_summary_tab), caption = caption)
}


#' S3 methods for printing model fit summary for glmfit class object
#'
#' Prints the modeling summary for the glm fit (\code{stats::glm.fit} or \code{speedglm::speedglm.wfit})
#' @param x The model fit object produced by functions stremr:::fit.glm or stremr:::fit.speedglm
#' @param ... Additional options passed on to \code{summary.GLMmodel}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary.GLMmodel}.
#' @export
print.GLMmodel <- function(x, ...) {
  model.summary <- summary(x, ...)
  cat(paste(model.summary, collapse = '\n'))
}

#' S3 methods for fit summary for glmfit class
#'
#' Prints the modeling summary for the GLM fit (\code{stats::glm.fit} or \code{speedglm::speedglm.wfit})
#' @param object The model fit object produced by functions stremr:::glmfit.glm or stremr:::glmfit.speedglm
#' @param format_table Format the coefficients into a data.frame table?
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.GLMmodel <- function(object, format_table = TRUE, ...) {
  makeModelCaption <- function(object) {
    return(
      "Model: " %+% object$params$outvar %+% " ~ " %+% paste0(object$params$predvars, collapse = " + ") %+% "; \\
       Stratify: " %+% object$params$stratify %+% "; \\
       N: " %+% prettyNum(object$nobs, big.mark = ",", scientific = FALSE) %+% "; \\
       Fit function: " %+% object$fitfunname
    )
  }
  nobs <- object$nobs
  coef_out <- object$coef
  if (format_table) {
    if (is.null(coef_out)) {
      coef_out <- "---"; names(coef_out) <- coef_out
    }
    coef_out <- data.table::data.table(Terms = names(coef_out), Coefficients = as.vector(coef_out))
    # coef_out <- data.frame(Terms = object$params$predvars, Coefficients = as.vector(coef_out))
    # rownames(coef_out) <- NULL
  }
  pander::set.caption(makeModelCaption(object))
  out <- pander::pander_return(coef_out, justify = c('right', 'left'))
  out
}

#' S3 methods for fit summary from xgboost
#'
#' Prints the modeling summary for the xgboost model fit (see \code{xgboost} R package).
#' @param object The model fit object produced by xgboost (and extracted with \code{getmodel_byname}).
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.xgb.Booster <- function(object, ...) {
  out <- NULL

    params <- object$params
    params <- lapply(params, function(arg) if (length(arg) > 1) {paste0(arg, collapse = ",")} else {arg})

    all_params <- t(data.table::as.data.table(params)) # all_params <- t(data.frame(params))
    all_params <- data.table::data.table(
      parameter = rownames(all_params),
      value = all_params[,1]) # , row.names = NULL

    all_params_pander <- pander::pander_return(all_params, caption  = "Model Parameters", justify = c('left', 'center'))
    out <- c(out, all_params_pander)

    training_stats <- data.table::data.table(
      name = c("best_iteration", "best_ntreelimit", "niter"),
      value = c(object$best_iteration, object$best_ntreelimit, object$niter)
      )

    training_stats_out <- pander::pander_return(training_stats, caption  = "Model Training Stats")
    out <- c(out, training_stats_out)

    # -----------------------------------------------------------------
    # model metrics (training and validation):
    # -----------------------------------------------------------------
    if (!is.null(object$best_score)) {
      metric_name <- attr(object$best_score, "names")
      performance <- data.table::data.table(
        metric_name = metric_name,
        best_score = object$best_score)

      performance_out <- pander::pander_return(performance, caption  = "Model Performance")
      out <- c(out, performance_out)
    }

    metrics <- object$evaluation_log
    if (nrow(metrics) > 10) metrics <- metrics[c(1:5, (nrow(metrics)-4):nrow(metrics)), ]
    model_metrics_out <- pander::pander_return(metrics, caption  = "Model Performance By Iteration")
    out <- c(out, model_metrics_out)

    # -----------------------------------------------------------------
    # variable importance:
    # -----------------------------------------------------------------
    feature_names <- object[["params"]][["feature_names"]]
    caption <- "Feature Importance"

    if (class(object) %in% "xgb.cv.synchronous") {
      object <- object[["models"]][[1]]
      caption <- caption %+% " for Training Fold 1"
    }

    try({
      importance_matrix <- xgboost::xgb.importance(feature_names = feature_names, model = object)
      importance_out <- pander::pander_return(importance_matrix, caption  = caption)
      out <- c(out, importance_out)
    }, silent = TRUE)

    ## To plot the feature importance:
    # xgboost::xgb.plot.importance(importance_matrix = importance_matrix)

  return(out)
}

#' @rdname summary.xgb.Booster
#' @export
summary.xgb.cv.synchronous <- function(object, ...) summary.xgb.Booster(object, ...)

#' S3 methods for fit summary for h2o
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param object The model fit object produced by h2o (and extracted with \code{getmodel_byname}).
#' @param only.coefs Skip all of the h2o-specific model stats, only print the coefficients table (when running \code{fit.algorithm = "glm"}).
# format_table Format the coefficients into a data.frame table (when running \code{fit.algorithm = "glm"})?
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.H2ORegressionModel <- function(object, only.coefs = FALSE, ...) {

  # object <- object$object.object
  modelID <- object@model$training_metrics@metrics$model$name
  out <- NULL

  # -----------------------------------------------------------------
  # some basic model info:
  # -----------------------------------------------------------------
  # coef_summary_out <- summary.GLMmodel(object, format_table)
  if (!is.null(object@model$coefficients_table)) {
    coef_summary_out <- pander::pander_return(object@model$coefficients_table, caption = attributes(object@model$coefficients_table)$description)
    out <- c(out, coef_summary_out)
  }

  if (!only.coefs) {
    # -----------------------------------------------------------------
    # model summary:
    # -----------------------------------------------------------------
    model_summary <- object@model$model_summary
    caption_summary <- attributes(model_summary)$header %+% " (Model ID: " %+% modelID %+%")"
    model_summary_out <- pander::pander_return(model_summary, caption = caption_summary)
    out <- c(out, model_summary_out)

    # -----------------------------------------------------------------
    # model parameters:
    # -----------------------------------------------------------------
    covars <- paste0(object@parameters$x, collapse = ",")
    predictors <- pander::pander_return(data.frame(predictors = covars))
    out <- c(out, predictors)

    params <- object@parameters[!names(object@parameters) %in% c("x", "model_id")]
    params <- lapply(params, function(arg) if (length(arg) > 1) {paste0(arg, collapse = ",")} else {arg})

    all_params <- t(data.table::as.data.table(params))

    all_params <- data.table::data.table(
      parameter = rownames(all_params),
      value = all_params[,1])

    # all_params <- data.frame(param = names(params), value = unlist(params), row.names = NULL)
    all_params_pander <- pander::pander_return(all_params, caption  = "Detailed Model Parameters", justify = c('left', 'center'))
    out <- c(out, all_params_pander)

    # -----------------------------------------------------------------
    # training data metrics:
    # -----------------------------------------------------------------
    train_model_metrics_out <- pander::pander_return(object@model$training_metrics, type = "Training")
    out <- c(out, train_model_metrics_out)

    # -----------------------------------------------------------------
    # validation data metrics:
    # -----------------------------------------------------------------
    H2OBinomialMetrics_val <- object@model$validation_metrics
      if (!is.null(H2OBinomialMetrics_val@metrics)) {
      valid_model_metrics_out <- pander::pander_return(H2OBinomialMetrics_val, type = "Validation")
      out <- c(out, valid_model_metrics_out)
    }

    # -----------------------------------------------------------------
    # cross validation data metrics:
    # -----------------------------------------------------------------
    H2OBinomialMetrics_xval <- object@model$cross_validation_metrics
    if (!is.null(H2OBinomialMetrics_xval@metrics)) {
      xval_model_metrics_out <- pander::pander_return(H2OBinomialMetrics_xval, type = "Cross-validation")
      out <- c(out, xval_model_metrics_out)
    }

    # -----------------------------------------------------------------
    # variable importance:
    # -----------------------------------------------------------------
    # h2o.varimp(object)
    var_imp <- object@model$variable_importances
    var_imp_cap <- attributes(var_imp)$header %+% "Model ID: " %+% modelID %+%")"
    var_imp_out <- pander::pander_return(var_imp, caption = var_imp_cap)
    out <- c(out, var_imp_out)
  }
  return(out)
}

#' @rdname summary.H2ORegressionModel
#' @export
summary.H2OBinomialModel <- function(object, only.coefs = FALSE, ...) summary.H2ORegressionModel(object, only.coefs, ...)


#' S3 methods for printing model fit summaries as pander tables
#'
#' Used internally to prints the modeling summary stats for reporting (see \code{\link{make_model_report}}).
#' @param model The model fit object produced by \code{get_fit} function.
#' @param only.coefs Set to \code{TRUE} to only print the table of coefficients (when using \code{fit.algorithm = "glm"}) and
#' omit all other model-related output statistics.
#' @param ... Additional options passed on to \code{summary(...)}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary(...)}.
#' @export
print_tables <- function(model, only.coefs = FALSE, ...) { UseMethod("print_tables") }

#' @rdname print_tables
#' @export
print_tables.face.sparse <- function(model, only.coefs = FALSE, ...) {
  cat("...face objects do not provide any model summary output...")
}

#' @rdname print_tables
#' @export
print_tables.brokenstick <- function(model, only.coefs = FALSE, ...) {
  old.opt <- pander::panderOptions('knitr.auto.asis')
  pander::panderOptions('knitr.auto.asis', TRUE)
  print(model)
  pander::panderOptions('knitr.auto.asis', old.opt)
  return(invisible(NULL))
}

#' @rdname print_tables
#' @export
print_tables.xgb.Booster <- function(model, only.coefs = FALSE, ...) {
  cat(paste(summary(model, only.coefs, ...), collapse = '\n'))
  # paste(summary(model, only.coefs, ...), collapse = '\n')
}
#' @rdname print_tables
#' @export
print_tables.xgb.cv.synchronous <- function(model, only.coefs = FALSE, ...) print_tables.xgb.Booster(model, only.coefs, ...)


#' @rdname print_tables
#' @export
print_tables.H2OBinomialModel <- function(model, only.coefs = FALSE, ...) {
  cat(paste(summary(model, only.coefs, ...), collapse = '\n'))
  # paste(summary(model, only.coefs, ...), collapse = '\n')
}

#' @rdname print_tables
#' @export
print_tables.H2ORegressionModel <- function(model, only.coefs = FALSE, ...) {
  cat(paste(summary(model, only.coefs, ...), collapse = '\n'))
  # paste(summary(model, only.coefs, ...), collapse = '\n')
}

#' @rdname print_tables
#' @export
print_tables.GLMmodel <- function(model, only.coefs = FALSE, ...) {
  cat(paste(summary(model, only.coefs, ...), collapse = '\n'))
}

CVmetrics_H2Obasemodel <- function(basemodelfit) {

  CV <- basemodelfit@model$cross_validation_metrics@metrics
  cap <- CV$description
  CV.metrics.tab <- data.frame(
              model = basemodelfit@model_id,
              CV$nobs,
              CV$MSE,
              CV$RMSE,
              CV$logloss,
              CV$r2,
              CV$AUC,
              CV$Gini,
              CV$mean_per_class_error,
              CV$model_category)

  return(CV.metrics.tab)
}

CVsummary_H2Obasemodel <- function(basemodelfit) {
  out <- NULL
  CVsummarytab <- basemodelfit@model$cross_validation_metrics_summary
  CVsummarytab <- CVsummarytab[, c("mean", "sd")]
  caption <- attributes(basemodelfit@model$cross_validation_metrics_summary)$header %+% ": " %+% basemodelfit@model_id
  out <- c(out, pander::pander_return(CVsummarytab, caption = caption))
  return(out)
}

#' S3 methods for getting model fit summary for H2Oensemblemodel class object
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param object The model fit object produced by h2oEnsemble package
#' @param only.coefs Skip all of the h2o-specific model stats, only print the coefficients table (when running \code{fit.algorithm = "glm"}).
#' @param format_table Format the coefficients into a data.frame table (when running \code{fit.algorithm = "glm"})?
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.H2Oensemblemodel <- function(object, only.coefs = FALSE, format_table = TRUE, ...) {
  # object <- object$H2O.model.object
  out <- NULL
  x <- object$x
  y <- object$y
  family <- object$family
  learner <- object$learner
  metalearner <- object$metalearner
  Vfolds <- object$cvControl$V
  seed <- object$seed

  # -----------------------------------------------------------------
  # some basic model info:
  # -----------------------------------------------------------------
  CV.descr <- data.frame(rbind(metalearner = metalearner, Vfolds = Vfolds, seed = seed, N.learners = length(learner)))
  names(CV.descr) <- c("")
  CV.descr.tab <- pander::pander_return(CV.descr, caption = "SuperLearner settings")
  out <- c(out, CV.descr.tab)

  coef_summary_out <- summary.GLMmodel(object, format_table)
  if (!is.null(object@model$coefficients_table)) {
    coef_summary_out <- pander::pander_return(object@model$coefficients_table, caption = attributes(object@model$coefficients_table)$description)
    out <- c(out, coef_summary_out)
  }

  out <- c(out, coef_summary_out)

  if (!only.coefs) {
    # -----------------------------------------------------------------
    # model summary:
    # -----------------------------------------------------------------
    # str(object$metafit@model$model_summary)
    model_summary <- object$metafit@model$model_summary
    caption_summary <- attributes(model_summary)$header %+% " Wrapper: " %+% metalearner
    model_summary_out <- pander::pander_return(model_summary, caption = caption_summary)
    out <- c(out, model_summary_out)

    # -----------------------------------------------------------------
    # CV metrics for each base learner
    # -----------------------------------------------------------------
    CVmetrics_tab <- NULL
    for (basemodel in object$basefits) {
      CVmetrics_tab <- rbind(CVmetrics_tab, CVmetrics_H2Obasemodel(basemodel))
      # out <- c(out, CVsummary_H2Obasemodel(basemodel))
    }
    cap <- basemodel@model$cross_validation_metrics@metrics$description
    out <- c(out, pander::pander_return(CVmetrics_tab, caption = cap))

  }
  return(out)
}

#' S3 methods for printing model fit summary for H2Omodel class object
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param x The model fit object produced by any stremr S3 function starting with \code{stremr:::H2Omodel.}
#' @param ... Additional options passed on to \code{summary.H2Omodel}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary.H2Omodel}.
#' @export
print.H2Oensemblemodel <- function(x, ...) {
  model.summary <- summary(x, ...)
  cat(paste(model.summary, collapse = '\n'))
}
