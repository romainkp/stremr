## ----------------------------------------------------------------------------------
## A trivial class for dealing with deterministic outcome modeling
## ----------------------------------------------------------------------------------
ModelDeterministic  <- R6Class(classname = "ModelDeterministic",
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
      self$gstar.Name <- reg$model_contrl[["gstar.Name"]]
      self$modelfit.g <- reg$model_contrl[["modelfit.g"]]
      self$intervened_type <- reg$model_contrl[["intervened_type"]]
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
      # self$n <- data$nobs
      # self$define.subset.idx(data)
      # private$probA1 <- data$get.outvar(TRUE, self$gstar.Name)
      # private$.isNA.probA1 <- is.na(private$probA1)
      # self$subset_idx <- rep.int(TRUE, self$n)
      # self$subset_idx <- seq_len(self$n)
      # private$.outvar <- data$get.outvar(TRUE, self$getoutvarnm) # Always a vector of 0/1
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

    predictAeqa = function(newdata, ...) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for Aobs[i] (n vector of a`s)
      assert_that(self$is.fitted)
      self$n <- newdata$nobs
      # self$define.subset.idx(newdata)
      if (missing(newdata)) {
        stop("newdata must be provided for evaluating gstar")
        # Aobs <- self$getoutvarval
      } else {
        Aobs <-  newdata$get.outvar(TRUE, self$getoutvarnm)
        Astar <- newdata$get.outvar(TRUE, self$gstar.Name)
      }

      # gstar <- rep.int(1L, self$n) # for missing values, the likelihood is always set to P(A = a) = 1.
      if (self$intervened_type %in% "bin") {
        # -------------------------------------------------
        # map gstar for intervention on binary A
        # -------------------------------------------------
        # check that observed exposure is always a vector of integers
        if(!is.integerish(Aobs)) {
          stop("Not possible to intervene on continuous node ('" %+% self$getoutvarnm %+% "') with a static intervention -- please change the intervention type by setting 'intervened_type_TRT' / 'intervened_type_MONITOR' to 'shift' or 'MSM'.")
        }
        gstar <- Astar^(Aobs) * (1 - Astar)^(1L - Aobs)
      } else if (self$intervened_type %in% "shift") {
        # -------------------------------------------------
        ## map gstar for a delta(W) shift of continuous A
        # -------------------------------------------------
        ## 1. define the new outcome node to evalute g(A-\delta(W)), defined as Anew = 2*Aobs - Astar
        newdata$dat.sVar[, "Anew.delta.shift.star.tmp" := 2*Aobs - Astar]
        ## 2. swap the names of gnode and self$gstar.Name
        newdata$swapNodes(current = self$getoutvarnm, target = "Anew.delta.shift.star.tmp")
        ## 3. call predict on original g fit, but using Anew.tmp as the outcome
        gstar <- self$modelfit.g$predictAeqa(newdata)
        ## 4. swap the node names back
        newdata$swapNodes(current = "Anew.delta.shift.star.tmp", target = self$getoutvarnm)
        newdata$dat.sVar[, "Anew.delta.shift.star.tmp" := NULL]
      } else if (self$intervened_type %in% "MSM") {
        # -------------------------------------------------
        ## map gstar for MSMs, which is just const 1
        # -------------------------------------------------
        gstar <- rep.int(1L, self$n)
      } else {
        stop("Unrecognized value for 'intervened_type_TRT' or 'intervened_type_MONITOR', these can only be 'bin', 'shift' or 'MSM', please consult the manual.")
      }

      self$wipe.alldat # to save RAM space when doing many stacked regressions wipe out all internal data:
      return(gstar)
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
