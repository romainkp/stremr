## ---------------------------------------------------------------------
#' R6 class for fitting and predicting joint probability for a univariate categorical summary A[j]
#'
#' This R6 class defines and fits a conditional probability model \code{P(A[j]|W,...)} for a univariate
#'  categorical summary measure \code{A[j]}. This class inherits from \code{\link{ModelGeneric}} class.
#'  Defines the fitting algorithm for a regression model \code{A[j] ~ W + ...}.
#'  Reconstructs the likelihood \code{P(A[j]=a[j]|W,...)} afterwards.
#'  Categorical \code{A[j]} is first redefined into \code{length(levels)} bin indicator variables, where
#'  \code{levels} is a numeric vector of all unique categories in \code{A[j]}.
#'  The fitting algorithm estimates the binary regressions for hazard for each bin indicator, \code{Bin_A[j][i] ~ W},
#'  i.e., the probability that categorical \code{A[j]} falls into bin \code{i}, \code{Bin_A[j]_i},
#'  given that \code{A[j]} does not fall in any prior bins \code{Bin_A[j]_1, ..., Bin_A[j]_{i-1}}.
#'  The dataset of bin indicators (\code{BinA[j]_1,...,BinA[j]_M}) is created
#'  inside the passed \code{data} or \code{newdata} object when defining \code{length(levels)} bins for \code{A[j]}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{levels}} - .
#' \item{\code{nbins}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, DataStorageClass.g0, ...)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
ModelCategorical <- R6Class(classname = "ModelCategorical",
  inherit = ModelGeneric,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the categorical outcome var (sA[j])
    levels = numeric(),       # all unique values for sA[j] sorted in increasing order
    nbins = integer(),
    bin_nms = character(),
    # Define settings for fitting cat sA and then call $new for super class (ModelGeneric)
    initialize = function(reg, DataStorageClass.g0, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      # Define the number of bins (no. of binary regressions to run) based on number of unique levels for categorical sVar:
      if (self$reg$get.reg$censoring & gvars$verbose) {
        message("...fitting a model for categorical censoring...")
      }
      if (gvars$verbose == 2) print("ModelCategorical outcome: "%+%self$outvar)

      assert_that(is.DataStorageClass(DataStorageClass.g0))
      self$levels <- DataStorageClass.g0$detect.cat.sVar.levels(reg$outvar)
      self$nbins <- length(self$levels)
      self$bin_nms <- DataStorageClass.g0$bin.nms.sVar(reg$outvar, self$nbins)

      # Instead of defining new RegressionClass, just clone the parent reg object and adjust the outcomes
      bin_regs <- self$reg$clone()
      bin_regs$outvar.class <- as.list(rep_len(gvars$sVartypes$bin, self$nbins))
      names(bin_regs$outvar.class) <- self$bin_nms
      bin_regs$outvar <- self$bin_nms
      bin_regs$predvars <- self$reg$predvars
      # subsetting variable names (always subset by non-missing values for the current bin column)
      bin_regs$subset_vars <- lapply(self$bin_nms, function(var) { c(var, self$reg$subset_vars)})
      names(bin_regs$subset_vars) <- self$bin_nms
      bin_regs$reg_hazard <- TRUE   # Don`t add degenerate bins as predictors in each binary regression

      super$initialize(reg = bin_regs, no_set_outvar = TRUE, ...)
    },

    # Transforms data for categorical outcome to bin indicators A[j] -> BinA[1], ..., BinA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in A - names have changed though)
    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      # Binirizes & saves binned matrix inside DataStorageClass for categorical sVar
      data$binirize.sVar(name.sVar = self$outvar, levels = self$levels)
      if (gvars$verbose == 2) {
        print("performing fitting for categorical outcome: " %+% self$outvar)
        print("freq counts by bin for categorical outcome: "); print(table(data$get.sVar(self$outvar)))
        print("binned dataset: "); print(head(cbind(sA = data$get.sVar(self$outvar), data$dat.bin.sVar), 5))
      }
      super$fit(data, ...) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")
      data$emptydat.bin.sVar # wiping out binirized mat in data object DataStorageClass...
      # self$wipe.alldat # wiping out all data traces in ContinModel...
      invisible(self)
    },

    # P(A=1|W=w): uses private$m.fit to generate predictions
    predict = function(newdata, ...) {
      # if (missing(newdata)) stop("must provide newdata")
      if (gvars$verbose == 2) print("performing prediction for categorical outcome: " %+% self$outvar)
      if (!missing(newdata)) assert_that(is.DataStorageClass(newdata))
      if (!missing(newdata)) newdata$binirize.sVar(name.sVar = self$outvar, levels = self$levels)
      super$predict(newdata, ...)
      if (!missing(newdata)) newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DataStorageClass object...
      invisible(self)
    },

    # Invisibly return cumm. prob P(A=a|W=w)
    # P(A=a|W=w) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    predictAeqa = function(newdata, ...) {
      if (gvars$verbose == 2) print("performing prediction for categorical outcome: " %+% self$outvar)
      if (!missing(newdata)) assert_that(is.DataStorageClass(newdata))
      if (!missing(newdata)) newdata$binirize.sVar(name.sVar = self$outvar, levels = self$levels)
      cumprodAeqa <- super$predictAeqa(newdata, ...)
      if (!missing(newdata)) newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      # self$wipe.alldat # wiping out all data traces in ModelCategorical
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)