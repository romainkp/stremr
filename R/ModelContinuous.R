ModelContinuous  <- R6Class(classname = "ModelContinuous",
  inherit = ModelCategorical,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    model_contrl = list(),
    models = list(),

    initialize = function(reg, ...) {
      model_contrl <- reg$model_contrl
      ## For continuous exposures, assume the user has already wrapped the sl3 binomial learner into correct condensier learners.
      if (!is.null(model_contrl[["models"]])) {
        self$models <- model_contrl[["models"]]
        model_contrl[["models"]] <- NULL
        if (!is(self$models, "Lrnr_base")) {
          stop("for continuous exposures, have to use sl3 package for defining learners / estimators")
        }
      } else {
        stop("for continuous exposures 'model' must be always always specified and it must be an sl3 learner object")
      }

      self$model_contrl <- model_contrl

      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar
      self$outvar.class <- reg$outvar.class
      self$outcome_type <- "continuous"

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
        print("New 'ModelContinuous' regression defined:"); print(self$show())
      }

      invisible(self)
    }
  )
)