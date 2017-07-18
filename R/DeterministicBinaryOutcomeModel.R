## ----------------------------------------------------------------------------------
## A trivial class for dealing with deterministic outcome modeling
## ----------------------------------------------------------------------------------
DeterministicBinaryOutcomeModel  <- R6Class(classname = "DeterministicBinaryOutcomeModel",
  inherit = BinaryOutcomeModel,
  cloneable = TRUE,
  portable = TRUE,
  class = TRUE,
  public = list(
    gstar.Name = character(),
    is.fitted = TRUE,

    initialize = function(reg, ...) {
      self$model_contrl <- reg$model_contrl
      self$gstar.Name <- reg$model_contrl[["gstar.Name"]]
      assert_that(!is.null(self$gstar.Name))
      assert_that(is.string(reg$outvar))
      self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)
      self$ReplMisVal0 <- reg$ReplMisVal0
      invisible(self)
    },

    # if (predict) then use the same data to make predictions for all obs in self$subset_idx;
    # store these predictions in private$probA1 and private$probAeqa
    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      self$define.subset.idx(data)
      private$probA1 <- data$get.outvar(TRUE, self$gstar.Name)
      # private$.isNA.probA1 <- is.na(private$probA1)
      # self$subset_idx <- rep.int(TRUE, self$n)
      self$subset_idx <- seq_len(self$n)
      private$.outvar <- data$get.outvar(TRUE, self$getoutvarnm) # Always a vector of 0/1
      # private$.isNA.outvar <- is.na(private$.outvar)
      self$is.fitted <- TRUE
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      # self$wipe.alldat
      # **********************************************************************
      invisible(self)
    },

    # get the fixed (known) the gstar P(A^*(t) = 1|W, bar{L(t)});
    # should be already saved earlier in private$probA1, so there is nothing to do here
    predict = function(newdata, ...) {
      assert_that(self$is.fitted)
      return(invisible(self))
    },

    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      if (missing(newdata)) {
        indA <- self$getoutvarval
      } else {
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      }
      assert_that(is.integerish(indA)) # check that observed exposure is always a vector of integers
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      # probA1 <- private$probA1[self$getsubset]
      probA1 <- private$probA1
      probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      self$wipe.alldat # to save RAM space when doing many stacked regressions wipe out all internal data:
      return(probAeqa)
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
      private$probA1 <- NULL
      private$probAeqa <- NULL
      private$.outvar <- NULL
      self$subset_idx <- NULL
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$subset_idx },
    getoutvarnm = function() { self$outvar },
    getoutvarval = function() { private$.outvar }
  ),

  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    .outvar = NULL,
    # .isNA.outvar = NULL,
    probA1 = NULL,      # Predicted probA^s=1 conditional on Xmat
    # .isNA.probA1 = NULL,
    probAeqa = NULL     # Likelihood of observing a particular value A^s=a^s conditional on Xmat
  )
)
