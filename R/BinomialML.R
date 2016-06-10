# IMPLEMENTING NEW CLASS FOR BINARY REGRESSION THAT DIRECTLY USES PRELOADED H20FRAME
# DOES NOT NEED TO RETRIEVE ANY COVARS, NEEDS TO ONLY EVALUATE THE SUBSET_IDX
# NEEDS TO be able to pass on THE REGRESSION SETTINGS FOR h2o-specific functions
# POSSIBLY NEEDS A SEPARATE RegressionClass
BinDatH2O  <- R6Class(classname = "BinDatH2O",
  inherit = BinDat,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(

    initialize = function(reg, ...) {
      assert_that(is.string(reg$outvar))
      assert_that(is.character(reg$predvars))
      self$outvar <- reg$outvar
      self$predvars <- reg$predvars

      self$subset_vars <- reg$subset_vars
      self$subset_expr <- reg$subset_exprs
      assert_that(length(self$subset_expr) <= 1)

      self$pool_cont <- reg$pool_cont
      self$outvars_to_pool <- reg$outvars_to_pool
      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- reg$nbins
      if (is.null(reg$subset_vars)) {self$subset_vars <- TRUE}
      assert_that(is.logical(self$subset_vars) || is.character(self$subset_vars)) # is.call(self$subset_vars) ||
      invisible(self)
    },

    # printing regression:
    show = function(print_format = TRUE) {
      if (print_format) {
        return("P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=", ") %+% ")" %+% ";\\ Stratify: " %+% self$subset_expr)
      } else {
        return(list(outvar = self$outvar, predvars = self$predvars, stratify = self$subset_expr))
      }
    },

    newdata = function(newdata, getoutvar = TRUE, ...) {
      assert_that(is.DataStorageClass(newdata))
      # CALL self$setdata.long() when: 1) self$pool_cont is TRUE & 2) more than one outvars_to_pool
      if (self$pool_cont && length(self$outvars_to_pool)>1) {
        self$setdata.long(data = newdata, ...)
      } else {
        self$setdata(data = newdata, getoutvar, ...)
      }
      invisible(self)
    },

    define.subset_idx = function(data) {
      if (is.logical(self$subset_vars)) {
        subset_idx <- self$subset_vars
      } else if (is.call(self$subset_vars)) {
        stop("calls aren't allowed in BinDat$subset_vars")
      } else if (is.character(self$subset_vars)) {
        subset_idx <- data$evalsubst(subset_vars = self$subset_vars, subset_expr = self$subset_expr)
      }
      assert_that(is.logical(subset_idx))
      if ((length(subset_idx) < self$n) && (length(subset_idx) > 1L)) {
        if (gvars$verbose) message("subset_idx has smaller length than self$n; repeating subset_idx p times, for p: " %+% data$p)
        subset_idx <- rep.int(subset_idx, data$p)
        if (length(subset_idx) != self$n) stop("BinDat$define.subset_idx: self$n is not equal to nobs*p!")
      }
      assert_that((length(subset_idx) == self$n) || (length(subset_idx) == 1L))
      return(subset_idx)
    },

    # Sets X_mat, Yvals, evaluates subset and performs correct subseting of data
    # everything is performed using data$ methods (data is of class DataStorageClass)
    setdata = function(data, getoutvar, ...) {
      assert_that(is.DataStorageClass(data))
      self$n <- data$nobs
      self$subset_idx <- self$define.subset_idx(data)
      if (getoutvar) private$Y_vals <- data$get.outvar(self$subset_idx, self$outvar) # Always a vector
      if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
        private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
        colnames(private$X_mat) <- c("Intercept", self$predvars)
      } else {
        # *** THIS IS THE ONLY LOCATION IN THE PACKAGE WHERE CALL TO DataStorageClass$get.dat.sVar() IS MADE ***
        if (length(self$predvars)==0L) {
          private$X_mat <- as.matrix(rep.int(1L, sum(self$subset_idx)), ncol=1)
        } else {
          private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.dat.sVar(self$subset_idx, self$predvars)))
        }
        colnames(private$X_mat)[1] <- "Intercept"
        # To find and replace misvals in X_mat:
        if (self$ReplMisVal0) private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace
      }

      self$DataStorageObject <- data

      invisible(self)
    },
    # Generic prediction fun for logistic regression coefs, predicts P(A = 1 | newXmat)
    # No need for S3 for now, until need different pred. funs for different classes
    # Does not handle cases with deterministic Anodes in the original data..
    logispredict = function(m.fit) {
      assert_that(!is.null(private$X_mat)); assert_that(!is.null(self$subset_idx))
      # Set to default missing value for A[i] degenerate/degerministic/misval:
      # Alternative, set to default replacement val: pAout <- rep.int(gvars$misXreplace, newdatsum_obj$n)
      pAout <- rep.int(gvars$misval, self$n)
      if (sum(self$subset_idx > 0)) {
        eta <- private$X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
        pAout[self$subset_idx] <- match.fun(FUN = m.fit$linkfun)(eta)
      }
      return(pAout)
    }
  ),
  active = list(
    wipe.alldat = function() {
      private$probA1 <- NULL
      private$probAeqa <- NULL
      self$bindat$emptydata
      self$bindat$emptyY
      self$bindat$emptySubset_idx
      self$bindat$emptyN
      return(self)
    },
    getfit = function() { private$m.fit },
    getprobA1 = function() { private$probA1 },
    getsubset = function() { self$bindat$subset_idx },
    getoutvarval = function() { self$bindat$getY },
    getoutvarnm = function() { self$bindat$outvar }
  ),
  private = list(
    m.fit = list(),   # the model fit (coefficients)
    probA1 = NULL,    # Predicted probA^s=1 conditional on X_mat
    probAeqa = NULL   # Likelihood of observing a particular value A^s=a^s conditional on X_mat
  )
)