
## ---------------------------------------------------------------------
# R6 class for fitting and predicting with several stratified models for a single outcome variable (conditional on some covariate values)
#
# This R6 class defines and fits a conditional probability model \code{P(A[j]|W,...)} for a summary \code{A[j]}.
# This class inherits from \code{\link{ModelGeneric}} class.
# The stratification criteria is determined by the \code{R} expression in the field \code{subset_exprs}.
#
ModelStratified <- R6Class(classname = "ModelStratified",
  inherit = ModelGeneric,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the binary outcome var (sA[j])
    # levels = numeric(),     # all unique values for sA[j] sorted in increasing order
    # nbins = integer(),
    subset_exprs = NULL,
    # Define settings for fitting cat sA and then call $new for super class (ModelGeneric)
    # We produce a regression class with 3 outvars (same) and 3 outvar.class (same)
    # The daughter regression clases resulting from this need to be of class StratifiedRegressionModelClass
    initialize = function(reg, DataStorageClass.g0, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$reg$subset_exprs) > 1L)
      assert_that(length(self$reg$outvar) == 1L)
      if (gvars$verbose == 2)  {
        print("ModelStratified outcome: "%+%self$outvar)
        print("ModelStratified expressions: ("%+% paste(self$subset_exprs, collapse=",") %+% ")")
      }
      # print("self$reg before stratification:"); self$reg$show()
      stratify_regs <- self$reg$clone()  # Instead of defining new RegressionClass, just clone the parent reg object and adjust the outcomes
      stratify_regs$outvar <- rep_len(self$reg$outvar, length(self$reg$subset_exprs))
      stratify_regs$outvar.class <- as.list(rep_len(self$reg$outvar.class, length(self$reg$subset_exprs)))
      names(stratify_regs$outvar.class) <- stratify_regs$outvar
      # by turning stratify_regs$subset_exprs into a list the subsetting will be performed on subset_exprs during next S3 dispatch on SummaryModel
      stratify_regs$subset_exprs <- as.list(stratify_regs$subset_exprs)
      names(stratify_regs$subset_exprs) <- stratify_regs$outvar
      stratify_regs$reg_hazard <- TRUE
      super$initialize(reg = stratify_regs, no_set_outvar = TRUE, DataStorageClass.g0 = DataStorageClass.g0, ...)
    },
    # Transforms data for categorical outcome to bin indicators sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      if (gvars$verbose == 2) {
        print("performing fitting for outcome based on stratified model for outcome: " %+% self$outvar)
        # print("following subsets are defined: "); print(table(data$get.sVar(self$outvar)))
      }
      super$fit(data, ...) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")
      invisible(self)
    },
    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata, ...) {
      if (gvars$verbose == 2) print("performing prediction for outcome based on stratified model: " %+% self$outvar)
      if (!missing(newdata)) assert_that(is.DataStorageClass(newdata))
      super$predict(newdata, ...)
      invisible(self)
    },
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    # P(A=a|W=w) - calculating the likelihood for obsdat.A[i] (n vector of a's):
    predictAeqa = function(newdata, ...) {
      if (gvars$verbose == 2) print("performing prediction for outcome based on stratified model: " %+% self$outvar)
      if (!missing(newdata)) assert_that(is.DataStorageClass(newdata))
      cumprodAeqa <- super$predictAeqa(newdata, ...)
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },
    predictgstar = function(newdata, ...) {
      if (gvars$verbose == 2) print("performing prediction for outcome based on stratified model: " %+% self$outvar)
      if (!missing(newdata)) assert_that(is.DataStorageClass(newdata))
      cumprodAeqa <- super$predictgstar(newdata, ...)
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    }
  ),
  active = list(
    # cats = function() {seq_len(self$reg$nbins)}
  )
)