#----------------------------------------------------------------------------------
# Classes that control modelling of the multivariate joint probability model P(sA|sW).
#----------------------------------------------------------------------------------

#' @importFrom assertthat assert_that
NULL

## ---------------------------------------------------------------------------------
## S3 constructors for the summary model classes:
## ---------------------------------------------------------------------------------
newsummarymodel <- function(regClass, reg, DataStorageClass.g0, ...) { UseMethod("newsummarymodel") }
## Summary model constructor for generic regression with multivariate outcome, but one set of predictors
newsummarymodel.generic <- function(regClass, reg, DataStorageClass.g0, ...) ModelGeneric$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
## Summary model constructor for stratification (by reg$subset_exprs):
newsummarymodel.stratify <- function(regClass, reg, DataStorageClass.g0, ...) ModelStratified$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
## Summary model constructor for constant outcome A (for now define the usual binomial model, in the future might return the outcome instead):
newsummarymodel.constant <- function(regClass, reg, ...) ModelBinomial$new(reg = reg, ...)
## Summary model constructor for binary outcome A:
newsummarymodel.binomial <- function(regClass, reg, ...) ModelBinomial$new(reg = reg, ...)
## Summary model constructor for categorical outcome A:
newsummarymodel.categorical <- function(regClass, reg, DataStorageClass.g0, ...) ModelCategorical$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
## Summary model constructor for continuous outcome A:
newsummarymodel.continuous <-  function(regClass, reg, DataStorageClass.g0, ...) ModelContinuous$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)

## Summary model constructor for Q-learning (sequential regression):
newsummarymodel.Qlearn <- function(regClass, reg, ...) ModelQlearn$new(reg = reg, ...)
newsummarymodel.SDRQlearn <- function(regClass, reg, ...) SDRModelQlearn$new(reg = reg, ...)
newsummarymodel.SplitCVSDRQlearn <- function(regClass, reg, ...) SplitCVSDRModelQlearn$new(reg = reg, ...)
newsummarymodel.SDRtransformQModel <- function(regClass, reg, ...) SDRtransformQModel$new(reg = reg, ...)

## A trivial class for dealing with NULL outcome modeling (when MONITOR and / or CENS aren't specified)
newsummarymodel.NULL <- function(regClass, reg, ...) ModelNULLOutcome$new(reg = reg, ...)

prettyprint_ModelGeneric <- function(self, reg, all.outvar.bin) {
  print("#----------------------------------------------------------------------------------")
  print("New instance of ModelGeneric:")
  print("#----------------------------------------------------------------------------------")
  # if ("ListOfRegressionForms" %in% class(reg$RegressionForms)) {
  if ("ListOfRegressionForms" %in% class(reg)) {
    # print("...ListOfRegressionForms..."); if (!is.null(names(reg$RegressionForms))) {print(names(reg$RegressionForms))}
    print("...ListOfRegressionForms..."); if (!is.null(names(reg))) {print(names(reg))}
    print("Outcomes: "); print(paste0(get_outvars(reg), collapse = ","))
  } else {
    print("Outcomes: " %+% paste(reg$outvar, collapse = ", "))
    print("Predictors: " %+% paste(reg$predvars, collapse = ", "))
    # print("Outcomes: " %+% paste(reg$RegressionForms$outvar, collapse = ", "))
    # print("Predictors: " %+% paste(reg$RegressionForms$predvars, collapse = ", "))
  }
  print("No. of regressions: " %+% self$n_regs)
  print("All outcomes binary? " %+% all.outvar.bin)
  print("#----------------------------------------------------------------------------------")
}

## ---------------------------------------------------------------------
# Generic R6 class for modeling (fitting and predicting) P(A=a|W=w) where A can be a multivariate (A[1], ..., A[k]) and each A[i] can be binary, categorical or continous
#
# This R6 class Class for defining, fitting and predicting the probability model
#  \code{P(A|W)} under \code{g_star} or under \code{g_0} for variables
#  (\code{A,W}). Defines and manages the factorization of the multivariate conditional
#  probability model \code{P(A=a|...)} into univariate regression models
#  \code{A[j] ~ A[j-1] + ... + A[1] + W}. The class \code{self$new} method automatically
#  figures out the correct joint probability factorization into univariate conditional
#  probabilities based on name ordering provided by (\code{A_nms}, \code{W_nms}).
#  When the outcome variable \code{A[j]} is binary, this class will automatically call
#  a new instance of \code{ModelBinomial} class.
#  Provide \code{self$fit()} function argument \code{data} as a \code{\link{DataStorageClass}} class object.
#  This data will be used for fitting the model \code{P(A|W)}.
#  Provide \code{self$fit()} function argument \code{newdata} (also as \code{DataStorageClass} class) for predictions of the type
#  \code{P(A=1|W=w)}, where \code{w} values are coming from \code{newdata} object.
#  Finally, provide \code{self$predictAeqa} function \code{newdata} argument
#  (also \code{DataStorageClass} class object) for getting the likelihood predictions \code{P(A=sa|W=w)}, where
#  both, \code{sa} and \code{sw} values are coming from \code{newdata} object.
#
# @docType class
# @format An \code{\link{R6Class}} generator object
# @keywords R6 class
# @details
# \itemize{
# \item{\code{n_regs}} - .
# }
# @section Methods:
# \describe{
#   \item{\code{new(reg, ...)}}{...}
#   \item{\code{length}}{...}
#   \item{\code{getPsAsW.models}}{...}
#   \item{\code{getcumprodAeqa}}{...}
#   \item{\code{copy.fit(ModelGeneric)}}{...}
#   \item{\code{fit(data)}}{...}
#   \item{\code{predict(newdata)}}{...}
#   \item{\code{predictAeqa(newdata, ...)}}{...}
# }
# @section Active Bindings:
# \describe{
#   \item{\code{wipe.alldat}}{...}
# }
ModelGeneric <- R6Class(classname = "ModelGeneric",
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),      # outcome name(s)
    predvars = character(),    # names of predictor vars
    n_regs = integer(),        # total no. of reg. models (logistic regressions)
    initialize = function(reg, no_set_outvar = FALSE, ...) {
      self$reg <- reg
      if (!no_set_outvar) self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      # Number of sep. regressions to run, based on the number of outcomes in SingleRegressionFormClass length(reg$outvar) or the number of objects in ListOfRegressionForms
      self$n_regs <- get_n_regs(reg)
      all.outvar.bin <- FALSE
      if (!("ListOfRegressionForms" %in% class(reg))) {
        all.outvar.bin <-  all(reg$outvar.class %in% gvars$sVartypes$bin)
      }
      if (gvars$verbose == 2) prettyprint_ModelGeneric(self, reg, all.outvar.bin)
      # Factorize the joint into univariate regressions, by dimensionality of the outcome variable (sA_nms):
      for (k_i in 1:self$n_regs) {

        if ("ListOfRegressionForms" %in% class(reg)) {
          reg_i <- reg[[k_i]]
        } else {
          reg_i <- reg$clone()
          reg_i$ChangeManyToOneRegresssion(k_i, reg)
        }
        # Calling the constructor for the summary model P(sA[j]|\bar{sA}[j-1], sW}), dispatching on reg_i class
        regS3class <- reg_i$S3class
        if (is.null(regS3class)) {
          regS3class <- "generic"; class(regS3class) <- "generic"
        }

        PsAsW.model <- newsummarymodel(regS3class, reg_i, ...)
        private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
        names(private$PsAsW.models)[k_i] <- "P(sA|sW)."%+%k_i
      }
      invisible(self)
    },
    length = function(){ base::length(private$PsAsW.models) },
    getPsAsW.models = function() { private$PsAsW.models },  # get all summary model objects (one model object per outcome var sA[j])
    # get joint prob as a vector of the cumulative prod over j for P(sA[j]=a[j]|W)
    getcumprodAeqa = function() { private$cumprodAeqa },  
    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      # serial loop over all regressions in PsAsW.models:
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$fit(data = data, ...)
      }
      invisible(self)
    },
    # P(A=1|W=w): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (!missing(newdata)) assert_that(is.DataStorageClass(newdata))
      # serial loop over all regressions in PsAsW.models:
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$predict(newdata = newdata)
      }
      invisible(self)
    },
    # WARNING: Next 2 methods cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Uses daughter objects (stored from prev call to fit()) to get joint predictions for P(A=obsdat.A|W=w) (factorized)
    # Invisibly returns the joint probability P(A=a|W=w), also aves it as a private field "cumprodAeqa"
    # P(A=a|W=w) - calculating the likelihood for obsdat.A[i] (n vector of a's):
    predictAeqa = function(newdata, n, ...) {
      if (!missing(newdata)) {
        assert_that(is.DataStorageClass(newdata))
        n <- newdata$nobs
      }
      jprodAeqa <- rep.int(1L, n)
      # loop over all regressions in PsAsW.models:
      for (k_i in seq_along(private$PsAsW.models)) {
        jprodAeqa <- jprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, n = n, ...)
      }
      private$cumprodAeqa <- jprodAeqa
      return(jprodAeqa)
    },
    # Same as above, but calling predictgstar(...) methods, plus not saving the joint prob
    predictgstar = function(newdata, n, ...) {
      if (!missing(newdata)) {
        assert_that(is.DataStorageClass(newdata))
        n <- newdata$nobs
      }
      jprodAeqa <- rep.int(1L, n)
      # loop over all regressions in PsAsW.models:
      for (k_i in seq_along(private$PsAsW.models)) {
        jprodAeqa <- jprodAeqa * private$PsAsW.models[[k_i]]$predictgstar(newdata = newdata, n = n, ...)
      }
      # private$jprodAeqa <- jprodAeqa
      return(jprodAeqa)
    },

    # get pre-saved predictions P(Q=1) from the K indexed model fit for unique n obs.
    # predictRegK = function(K, n) {
    #   if (length(private$PsAsW.models) < K) stop("invalid arg K; this object contains only " %+% length(private$PsAsW.models) %+% " different model fits.")
    #   return(private$PsAsW.models[[K]]$predictAeqa(n = n))
    # },

    # call itself until reaches a terminal model fit with coefficients + regression returned with show()
    get.fits = function() {
      res_models <- NULL
      for (k_i in seq_along(private$PsAsW.models)) {
        res <- private$PsAsW.models[[k_i]]$get.fits()
        if (is.list(res)) res_models <- c(res_models, res)
      }
      return(res_models)
    },
    # call itself until reaches a terminal model fit with coefficients + regression returned with show()
    get.model.summaries = function() {
      res_models <- NULL
      for (k_i in seq_along(private$PsAsW.models)) {
        res <- private$PsAsW.models[[k_i]]$get.model.summaries()
        if (is.list(res)) res_models <- c(res_models, res)
      }
      return(res_models)
    }

  ),
  active = list(
    # recursively call all saved daughter model fits and wipe out any traces of saved data
    wipe.alldat = function() {
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$wipe.alldat
      }
      return(self)
    }
  ),
  private = list(
    deep_clone = function(name, value) {
      # if value is is an environment, quick way to copy:
      # list2env(as.list.environment(value, all.names = TRUE), parent = emptyenv())
      # if a list of R6 objects, make a deep copy of each:
      if (name == "PsAsW.models") {
        lapply(value, function(PsAsW.model) PsAsW.model$clone(deep=TRUE))
      } else if (inherits(value, "R6")) { # to check the value is an R6 object:
        value$clone(deep=TRUE)
      } else {
        value # For all other fields, just return the value
      }
    },
    PsAsW.models = list(),
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)