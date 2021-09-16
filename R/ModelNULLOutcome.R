## ----------------------------------------------------------------------------------
## A trivial class for dealing with NULL outcome modeling (when MONITOR and / or CENS aren't specified)
## This class does nothing but simply returns a vector of (1,1,1,...) when predict methods are called.
## ----------------------------------------------------------------------------------
ModelNULLOutcome  <- R6Class(classname = "ModelNULLOutcome",
  inherit = ModelBinomial,
  cloneable = TRUE,
  portable = TRUE,
  class = TRUE,
  public = list(
    gstar.Name = character(),
    modelfit.g = NULL,
    intervened_type = "bin",
    is.fitted = TRUE,
    initialize = function(reg, ...) {
      self$model_contrl <- reg$model_contrl
      assert_that(is.null(reg$outvar) || reg$outvar == "NULL")
      self$outvar <- reg$outvar
      self$subset_vars <- reg$subset_vars
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$subset_exprs) <= 1)
      self$ReplMisVal0 <- reg$ReplMisVal0
      invisible(self)
    },

    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      self$n <- data$nobs
      invisible(self)
    },

    predict = function(newdata, ...) {
      return(invisible(self))
    },

    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      return(probAeqa)
    },

    predictgstar = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      probAeqa <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      return(probAeqa)
    },

    # Output info on the general type of regression being fitted:
    show = function(print_format = TRUE) {
      if (print_format) {
        return("P(" %+% self$outvar %+% "|" %+% "..." %+% ")")
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_exprs))
      }
    }
  ),

  active = list(
    wipe.alldat = function() {
      return(self)
    }
  )
)