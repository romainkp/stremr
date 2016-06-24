#----------------------------------------------------------------------------------
# Classes for modelling regression models with binary outcome Bin ~ Xmat
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

#' S3 methods for printing model fit summary for glmfit class object
#'
#' Prints the modeling summary for the GLM fit (\code{stats::glm.fit} or \code{speedglm::speedglm.wfit})
#' @param model.fit The model fit object produced by functions stremr:::glmfit.glm or stremr:::glmfit.speedglm
#' @param ... Additional options passed on to \code{summary.glmfit}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary.glmfit}.
#' @export
print.glmfit <- function(model.fit, ...) {
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
summary.glmfit <- function(model.fit, format_table = TRUE, ...) {
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

#' S3 methods for getting model fit summary for h2ofit class object
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param model.fit The model fit object produced by any stremr S3 function starting with \code{stremr:::h2ofit.}
#' @param only.coefs Skip all of the h2o-specific model stats, only print the coefficients table (when running \code{fit.algorithm = "GLM"}).
#' @param format_table Format the coefficients into a data.frame table (when running \code{fit.algorithm = "GLM"})?
#' @param ... Additional options (not used)
#' @return The markdown-formated model summary returned by \code{pander::pander_return}.
#' @export
summary.h2ofit <- function(model.fit, only.coefs = FALSE, format_table = TRUE, ...) {
  h2o.model <- model.fit$H2O.model.object
  modelID <- h2o.model@model$training_metrics@metrics$model$name
  out <- NULL

  # -----------------------------------------------------------------
  # some basic model info:
  # -----------------------------------------------------------------
  coef_summary_out <- summary.glmfit(model.fit, format_table)
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

#' S3 methods for printing model fit summary for h2ofit class object
#'
#' Prints the modeling summary for the h2o model fit (see \code{h2o} R package).
#' @param model.fit The model fit object produced by any stremr S3 function starting with \code{stremr:::h2ofit.}
#' @param only.coefs Skip all of the h2o-specific model stats, only print the coefficients table (when running \code{fit.algorithm = "GLM"}).
#' @param ... Additional options passed on to \code{summary.h2ofit}.
#' @return The output is printed with \code{cat}. To capture the markdown-formated model summary use \code{summary.h2ofit}.
#' @export
print.h2ofit <- function(model.fit, only.coefs = FALSE, ...) {
  model.summary <- summary(model.fit, only.coefs, ...)
  cat(paste(model.summary, collapse = '\n'))
}

## ---------------------------------------------------------------------
#' R6 class for fitting and making predictions for a single binary outcome regression model P(B | PredVars)
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
#' \item{cont.sVar.flag} - Is the original outcome variable continuous?
#' \item{bw.j} - Bin width (interval length) for an outcome that is a bin indicator of a discretized continous outcome.
#' \item{GLMpackage} - Controls which package will be used for performing model fits (\code{glm} or \code{speedglm}).
#' \item{binomialModelObj} - Pointer to an instance of \code{binomialModelObj} class that contains the data.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg)}}{Uses \code{reg} R6 \code{\link{RegressionClass}} object to instantiate a new model for a
#'   logistic regression with binary outcome.}
#'   \item{\code{show()}}{Print information on outcome and predictor names used in this regression model}
#'   \item{\code{fit()}}{...}
#'   \item{\code{copy.fit()}}{...}
#'   \item{\code{predict()}}{...}
#'   \item{\code{copy.predict()}}{...}
#'   \item{\code{predictAeqa()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{getoutvarnm}}{...}
#'   \item{\code{getoutvarval}}{...}
#'   \item{\code{getsubset}}{...}
#'   \item{\code{getprobA1}}{...}
#'   \item{\code{getfit}}{...}
#'   \item{\code{wipe.alldat}}{...}
#' }
#' @importFrom assertthat assert_that is.flag
#' @export
BinaryOutcomeModel  <- R6Class(classname = "BinaryOutcomeModel",
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    cont.sVar.flag = logical(),
    bw.j = numeric(),
    is.fitted = FALSE,

    binomialModelObj = NULL, # object of class binomialModelObj that is used in fitting / prediction, never saved (need to be initialized with $new())
    fit.package = c("speedglm", "glm", "h2o"),
    fit.algorithm = c("GLM", "GBM", "RF", "SL"),
    model_contrl = list(),

    n = NA_integer_,        # number of rows in the input data
    nbins = integer(),
    subset_vars = NULL,     # THE VAR NAMES WHICH WILL BE TESTED FOR MISSINGNESS AND WILL DEFINE SUBSETTING
    subset_exprs = NULL,     # THE LOGICAL EXPRESSION (ONE) TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)

    ReplMisVal0 = logical(),

    initialize = function(reg, ...) {
      model_contrl <- reg$model_contrl
      self$model_contrl <- reg$model_contrl
      if ("fit.package" %in% names(model_contrl)) {
        self$fit.package <- model_contrl[['fit.package']]
        assert_that(is.character(self$fit.package))
        assert_that(self$fit.package %in% c("speedglm", "glm", "h2o"))
      } else {
        self$fit.package <- reg$fit.package[1]
      }

      if ("fit.algorithm" %in% names(model_contrl)) {
        self$fit.algorithm <- model_contrl[['fit.algorithm']]
        assert_that(is.character(self$fit.algorithm))
        assert_that(self$fit.algorithm %in% c("GLM", "GBM", "RF", "SL"))
      } else {
        self$fit.algorithm <- reg$fit.algorithm[1]
      }

      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar

      assert_that(is.character(reg$predvars))
      self$predvars <- reg$predvars

      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)

      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- reg$nbins

      if (is.null(reg$subset_vars)) {self$subset_vars <- TRUE}
      assert_that(is.logical(self$subset_vars) || is.character(self$subset_vars)) # is.call(self$subset_vars) ||

      # ***************************************************************************
      # Add any additional options passed on to modeling functions as extra args
      # ***************************************************************************
      if (self$fit.package %in% "h2o") {
        self$binomialModelObj <- BinomialH2O$new(fit.algorithm = self$fit.algorithm, fit.package = self$fit.package, ParentModel = self, ...)
      } else {
        self$binomialModelObj <- BinomialGLM$new(fit.algorithm = self$fit.algorithm, fit.package = self$fit.package, ParentModel = self, ...)
      }

      if (gvars$verbose) {
        print("New instance of " %+% class(self)[1] %+% " :"); print(self$show())
      }
      # Get the bin width (interval length) for the current bin name self$getoutvarnm (for discretized continuous sA only):
      self$cont.sVar.flag <- self$getoutvarnm %in% names(reg$intrvls.width)
      if (self$cont.sVar.flag) {
        intrvl.idx <- which(names(reg$intrvls.width) %in% self$getoutvarnm)
        if (length(intrvl.idx) > 1) stop("non-unique names for intrvls.width in RegressionClass")
        self$bw.j <- reg$intrvls.width[intrvl.idx]
      } else {
        self$bw.j <- 1L
      }
      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and private$probAeqa
    fit = function(overwrite = FALSE, data, predict = FALSE, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      if (gvars$verbose) print("fitting the model: " %+% self$show())
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked

      self$define.subset.idx(data)

      private$model.fit <- self$binomialModelObj$fit(data, self$outvar, self$predvars, self$subset_idx, ...)
      private$model.fit$params <- self$show(print_format = FALSE)

      self$is.fitted <- TRUE
      if (predict) {
        self$predictAeqa(...)
      }

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$model.fit to generate predictions for data:
    predict = function(newdata, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata) && is.null(private$probA1)) {
        # stop("must provide newdata for BinaryOutcomeModel$predict()")
        private$probA1 <- self$binomialModelObj$predictP1(subset_idx = self$subset_idx)
      } else {
        self$n <- newdata$nobs
        self$define.subset.idx(newdata)
        private$probA1 <- self$binomialModelObj$predictP1(data = newdata, subset_idx = self$subset_idx)
      }
      return(invisible(self))
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, bw.j.sA_diff, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      }

      self$predict(newdata)

      if (missing(newdata)) {
        indA <- self$getoutvarval
      } else {
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      }

      assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      probA1 <- private$probA1[self$getsubset]
      assert_that(!any(is.na(probA1))) # check that predictions P(A=1 | dmat) exist for all obs.
      # Discrete version for the joint density:
      probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      # continuous version for the joint density:
      # probAeqa[self$getsubset] <- (probA1^indA) * exp(-probA1)^(1 - indA)
      # Alternative intergrating the last hazard chunk up to x:
      # difference of sA value and its left most bin cutoff: x - b_{j-1}
      if (!missing(bw.j.sA_diff)) {
        # + integrating the constant hazard all the way up to value of each sa:
        # probAeqa[self$getsubset] <- probAeqa[self$getsubset] * (1 - bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
        # cont. version of above:
        probAeqa[self$getsubset] <- probAeqa[self$getsubset] * exp(-bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
      }
      private$probAeqa <- probAeqa

      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      return(probAeqa)
    },

    sampleA = function(newdata, bw.j.sA_diff) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      assert_that(!missing(newdata))
      # Don't want to subset by the outvar, since binarized mat for cat outcome is not re-created when just sampling
      # But need to reset it back when done
      temp_subset_vars <- self$subset_vars
      self$subset_vars <- self$subset_vars[!self$subset_vars %in% self$outvar]

      self$predict(newdata)

      sampleA <- rep.int(0L, n)
      sampleA[self$getsubset] <- rbinom(n = n, size = 1, prob = probA1)
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      self$subset_vars <- temp_subset_vars
      # **********************************************************************
      return(sampleA)
    },

    define.subset.idx = function(data) {
      if (is.logical(self$subset_vars)) {
        subset_idx <- self$subset_vars
      } else if (is.call(self$subset_vars)) {
        stop("calls aren't allowed in binomialModelObj$subset_vars")
      } else if (is.character(self$subset_vars)) {
        subset_idx <- data$evalsubst(subset_vars = self$subset_vars, subset_exprs = self$subset_exprs)
      }
      assert_that(is.logical(subset_idx))
      if ((length(subset_idx) < self$n) && (length(subset_idx) > 1L)) {
        if (gvars$verbose) message("subset_idx has smaller length than self$n; repeating subset_idx p times, for p: " %+% data$p)
        subset_idx <- rep.int(subset_idx, data$p)
        if (length(subset_idx) != self$n) stop("binomialModelObj$define.subset.idx: self$n is not equal to nobs*p!")
      }
      assert_that((length(subset_idx) == self$n) || (length(subset_idx) == 1L))
      self$subset_idx <- subset_idx
      return(invisible(self))
    },

    # take fitted BinaryOutcomeModel class object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("BinaryOutcomeModel" %in% class(bin.out.model))
      private$model.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # take BinaryOutcomeModel class object that contains the predictions for P(A=1|sW) and save these predictions to self$
    copy.predict = function(bin.out.model) {
      assert_that("BinaryOutcomeModel" %in% class(bin.out.model))
      assert_that(self$is.fitted)
      private$probA1 <- bin.out.model$getprobA1
    },

    # Returns the object that contains the actual model fits (itself)
    get.fits = function() {
      model.fit <- self$getfit
      return(list(model.fit))
    },

    # Output info on the general type of regression being fitted:
    show = function(print_format = TRUE) {
      if (print_format) {
        return("P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=", ") %+% ")" %+% ";\\ Stratify: " %+% self$subset_exprs)
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs))
      }
    }
  ),

  active = list(
    wipe.alldat = function() {
      # private$probA1 <- NULL
      # private$probAeqa <- NULL
      self$subset_idx <- NULL
      self$binomialModelObj$emptydata
      self$binomialModelObj$emptyY
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    getoutvarval = function() { self$binomialModelObj$getY }
  ),
  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    probA1 = NULL,    # Predicted probA^s=1 conditional on Xmat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)
