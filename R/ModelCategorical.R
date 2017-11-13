# ## ---------------------------------------------------------------------
# #' R6 class for fitting and predicting joint probability for categorical outcome
# #'
# #' This R6 class defines and fits a conditional probability model \code{P(A[j]|W,...)} for 
# #'  categorical summary measure \code{A[j]}. This class inherits from \code{\link{ModelGeneric}} class.
# #'  Defines the fitting algorithm for a regression model \code{A[j] ~ W + ...}.
# #'  Reconstructs the likelihood \code{P(A[j]=a[j]|W,...)} afterwards.
# #'  Categorical \code{A[j]} is first redefined into \code{length(levels)} bin indicator variables, where
# #'  \code{levels} is a numeric vector of all unique categories in \code{A[j]}.
# #'  The fitting algorithm estimates the binary regressions for hazard for each bin indicator, \code{Bin_A[j][i] ~ W},
# #'  i.e., the probability that categorical \code{A[j]} falls into bin \code{i}, \code{Bin_A[j]_i},
# #'  given that \code{A[j]} does not fall in any prior bins \code{Bin_A[j]_1, ..., Bin_A[j]_{i-1}}.
# #'  The dataset of bin indicators (\code{BinA[j]_1,...,BinA[j]_M}) is created
# #'  inside the passed \code{data} or \code{newdata} object when defining \code{length(levels)} bins for \code{A[j]}.
# #'
# #' @docType class
# #' @format An \code{\link{R6Class}} generator object
# #' @keywords R6 class
# #' @details
# #' \itemize{
# #' \item{\code{reg}} - .
# #' \item{\code{outvar}} - .
# #' \item{\code{levels}} - .
# #' \item{\code{nbins}} - .
# #' }
# #' @section Methods:
# #' \describe{
# #'   \item{\code{new(reg, DataStorageClass.g0, ...)}}{...}
# #'   \item{\code{fit(data)}}{...}
# #'   \item{\code{predict(newdata)}}{...}
# #'   \item{\code{predictAeqa(newdata)}}{...}
# #' }
# #' @section Active Bindings:
# #' \describe{
# #'   \item{\code{cats}}{...}
# #' }
# #' @export
ModelCategorical  <- R6Class(classname = "ModelCategorical",
  inherit = ModelBinomial,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    model_contrl = list(),
    models = list(),

    initialize = function(reg, ...) {
      model_contrl <- reg$model_contrl

      if (!is.null(model_contrl[["models"]])) {
        self$models <- model_contrl[["models"]]
        model_contrl[["models"]] <- NULL
        if (!is(self$models, "Lrnr_base")) {
          stop("for categorical exposures, have to use sl3 package for defining learners / estimators")
        }

        ## For categorical exposure we always wrap the sl3 learner into condensier learner.
        ## This will factorize the likelihood appropriately, 
        ## then fit a separate logistic regression problem to each level of the categorical exposure, 
        ## and will perform the final predictions in the rigth way (as likelihood).
        self$models <- sl3::Lrnr_condensier$new(bin_estimator = self$models)

      } else {
        stop("for categorical exposures 'model' must be always always specified and it must be an sl3 learner object")
      }

      self$model_contrl <- model_contrl

      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar
      self$outvar.class <- reg$outvar.class

      assert_that(is.character(reg$predvars))
      self$predvars <- reg$predvars

      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)

      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- reg$nbins

      if (is.null(reg$subset_vars)) {self$subset_vars <- TRUE}
      assert_that(is.logical(self$subset_vars) || is.character(self$subset_vars)) # is.call(self$subset_vars) ||

      if (gvars$verbose) {
        print("New 'ModelCategorical' regression defined:"); print(self$show())
      }

      invisible(self)
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, indA, holdout = FALSE, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      if (missing(newdata) && !is.null(private$probAeqa)) {
        return(private$probAeqa)
      }
      self$predict(newdata, holdout)
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      probA1 <- private$probA1 # [self$getsubset]
      # check that predictions P(A=1 | dmat) exist for all obs (not NA)
      if (any(is.na(probA1))) stop("some of the modeling predictions resulted in NAs, which indicates an error in a prediction routine")
      ## Predictions for categorical / continuous are already reported as likelihood P(A=a|W) when using condensier
      if (is.list(probA1) || is.data.table(probA1) || is.data.frame(probA1)) {
        probA1 <- probA1[[1]]
      }
      likelihood <- probA1
      probAeqa[self$getsubset] <- likelihood
      private$probAeqa <- probAeqa
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # **********************************************************************
      return(probAeqa)
    }
  )
)