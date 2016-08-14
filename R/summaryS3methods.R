#----------------------------------------------------------------------------------
# S3 classes for printing model summaries
#----------------------------------------------------------------------------------

#' @importFrom pander pander
NULL

#' Pander method for H2OBinomialMetrics class
#'
#' Prints a H2OBinomialMetrics object in Pandoc's markdown.
#' @param H2OBinomialMetricsObject H2OBinomialMetrics object
#' @return By default this function outputs (see: \code{?cat}) the result.
#' If you would want to catch the result instead, then call \code{pander_return} instead.
#' @export
pander.H2OBinomialMetrics <- function(H2OBinomialMetricsObject) {
  modelID <- H2OBinomialMetricsObject@metrics$model$name
  metricsDF <- data.frame(
      metric = c(
          'MSE:       ',
          'R^2:       ',
          "LogLoss:   ",
          "AUC:       ",
          "Gini:      "),
      value = c(
          H2OBinomialMetricsObject@metrics$MSE,
          H2OBinomialMetricsObject@metrics$r2,
          H2OBinomialMetricsObject@metrics$logloss,
          H2OBinomialMetricsObject@metrics$AUC,
          H2OBinomialMetricsObject@metrics$Gini),
      stringsAsFactors = FALSE)

    if (H2OBinomialMetricsObject@algorithm == "glm") {
      glmDF <- data.frame(
        metric =
          c("Null Deviance:      ",
            "Residual Deviance:  ",
            "AIC:                "),
        value = c(H2OBinomialMetricsObject@metrics$null_deviance,
              H2OBinomialMetricsObject@metrics$residual_deviance,
              H2OBinomialMetricsObject@metrics$AIC),
        stringsAsFactors = FALSE)
      metricsDF <- rbind(metricsDF, glmDF)
    }
    colnames(metricsDF) <- NULL
    pander::pander(metricsDF, justify = c('left', 'center'), caption = "Model ID: " %+% modelID)

    cm <- h2o::h2o.confusionMatrix(H2OBinomialMetricsObject)
    if( !is.null(cm) ) {
      pander::pander(cm, caption = "Confusion Matrix for F1-optimal threshold" %+% " (Model ID: " %+% modelID %+%")" )
    }
    max_matrics <- H2OBinomialMetricsObject@metrics$max_criteria_and_metric_scores
    caption <- attributes(max_matrics)$header %+% ": " %+% attributes(max_matrics)$description
    pander::pander(max_matrics, caption = caption %+% " (Model ID: " %+% modelID %+%")")
    # attributes(max_matrics)$formats
    # cat("\nGains/Lift Table: Extract with `h2o.gainsLift(<model>, <data>)` or `h2o.gainsLift(<model>, valid=<T/F>, xval=<T/F>)`")
    if (!is.null(H2OBinomialMetricsObject@metrics$gains_lift_table)) {
     # print(H2OBinomialMetricsObject@metrics$gains_lift_table)
     gain_tab <- H2OBinomialMetricsObject@metrics$gains_lift_table
     pander::pander(gain_tab, caption = attributes(gain_tab)$header %+% " (Model ID: " %+% modelID %+%")")
    }

  return(invisible(H2OBinomialMetricsObject))
}

#' Pander method for H2ORegressionMetrics class
#'
#' Prints a H2ORegressionMetrics object in Pandoc's markdown.
#' @param H2ORegressionMetricsObject H2ORegressionMetrics object
#' @return By default this function outputs (see: \code{?cat}) the result.
#' If you would want to catch the result instead, then call \code{pander_return} instead.
#' @export
pander.H2ORegressionMetrics <- function(H2ORegressionMetricsObject) {
  return(NULL)
}

#' S3 methods for printing model fit summary for glmfit class object
#'
#' Prints the modeling summary for the glm fit (\code{stats::glm.fit} or \code{speedglm::speedglm.wfit})
#' @param model.fit The model fit object produced by functions stremr:::fit.glm or stremr:::fit.speedglm
#' @param ... Additional options passed on to \code{summary.GLMmodel}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary.GLMmodel}.
#' @export
print.GLMmodel <- function(model.fit, ...) {
  model.summary <- summary(model.fit, ...)
  cat(paste(model.summary, collapse = '\n'))
}

#' S3 methods for getting model fit summary for glmfit class object
#'
#' Prints the modeling summary for the GLM fit (\code{stats::glm.fit} or \code{speedglm::speedglm.wfit})
#' @param model.fit The model fit object produced by functions stremr:::glmfit.glm or stremr:::glmfit.speedglm
#' @param format_table Format the coefficients into a data.frame table?
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.GLMmodel <- function(model.fit, format_table = TRUE, ...) {
  makeModelCaption <- function(model.fit) {
    return(
      "Model: " %+% model.fit$params$outvar %+% " ~ " %+% paste0(model.fit$params$predvars, collapse = " + ") %+% "; \\
       Stratify: " %+% model.fit$params$stratify %+% "; \\
       N: " %+% prettyNum(model.fit$nobs, big.mark = ",", scientific = FALSE) %+% "; \\
       Fit function: " %+% model.fit$fitfunname
    )
  }
  nobs <- model.fit$nobs
  coef_out <- model.fit$coef
  if (format_table) {
    if (is.null(coef_out)) {
      coef_out <- "---"; names(coef_out) <- coef_out
    }
    coef_out <- data.frame(Terms = names(coef_out), Coefficients = as.vector(coef_out))
    rownames(coef_out) <- NULL
  }
  pander::set.caption(makeModelCaption(model.fit))
  out <- pander::pander_return(coef_out, justify = c('right', 'left'))
  out
}

#' S3 methods for getting model fit summary for H2Omodel class object
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param model.fit The model fit object produced by any stremr S3 function based on h2o
#' @param only.coefs Skip all of the h2o-specific model stats, only print the coefficients table (when running \code{fit.algorithm = "glm"}).
#' @param format_table Format the coefficients into a data.frame table (when running \code{fit.algorithm = "glm"})?
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.H2Omodel <- function(model.fit, only.coefs = FALSE, format_table = TRUE, ...) {
  h2o.model <- model.fit$H2O.model.object
  modelID <- h2o.model@model$training_metrics@metrics$model$name
  out <- NULL

  # -----------------------------------------------------------------
  # some basic model info:
  # -----------------------------------------------------------------
  coef_summary_out <- summary.GLMmodel(model.fit, format_table)
  out <- c(out, coef_summary_out)

  if (!only.coefs) {
    # -----------------------------------------------------------------
    # model summary:
    # -----------------------------------------------------------------
    model_summary <- h2o.model@model$model_summary
    caption_summary <- attributes(model_summary)$header %+% " (Model ID: " %+% modelID %+%")"
    model_summary_out <- pander::pander_return(model_summary, caption = caption_summary)
    out <- c(out, model_summary_out)

    # -----------------------------------------------------------------
    # training data metrics:
    # -----------------------------------------------------------------
    H2OBinomialMetrics_training <- h2o.model@model$training_metrics
    train_model_metrics_out <- pander::pander_return(H2OBinomialMetrics_training)
    out <- c(out, train_model_metrics_out)

    # -----------------------------------------------------------------
    # variable importance:
    # -----------------------------------------------------------------
    var_imp <- h2o.model@model$variable_importances
    var_imp_cap <- attributes(var_imp)$header %+% "Model ID: " %+% modelID %+%")"
    var_imp_out <- pander::pander_return(var_imp, caption = var_imp_cap)
    out <- c(out, var_imp_out)
  }
  return(out)
}

#' S3 methods for printing model fit summary for H2Omodel class object
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param model.fit The model fit object produced by any stremr S3 function starting with \code{stremr:::H2Omodel.}
#' @param only.coefs Skip all of the h2o-specific model stats, only print the coefficients table (when running \code{fit.algorithm = "GLM"}).
#' @param ... Additional options passed on to \code{summary.H2Omodel}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary.H2Omodel}.
#' @export
print.H2Omodel <- function(model.fit, only.coefs = FALSE, ...) {
  model.summary <- summary(model.fit, only.coefs, ...)
  cat(paste(model.summary, collapse = '\n'))
}