#' @import data.table
NULL

# Generic for fitting the logistic (binomial family) GLM model
logisfit <- function(BinDatObject, ...) UseMethod("logisfit")

# S3 method for glm binomial family fit, takes BinDat data object:
logisfit.glmS3 <- function(BinDatObject, ...) {
  if (gvars$verbose) print("calling glm.fit...")
  Xmat <- BinDatObject$getXmat
  Y_vals <- BinDatObject$getY
  # Xmat has 0 rows: return NA's and avoid throwing exception:
  if (nrow(Xmat) == 0L) {
    m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    ctrl <- glm.control(trace = FALSE)
    SuppressGivenWarnings({
      m.fit <- stats::glm.fit(x = Xmat, y = Y_vals, family = binomial() , control = ctrl)
    }, GetWarningsToSuppress())
  }
  fit <- list(coef = m.fit$coef, linkfun = "logit_linkinv", fitfunname = "glm", nobs = nrow(Xmat))
  if (gvars$verbose) print(fit$coef)
  class(fit) <- c(class(fit), c("glmS3"))
  return(fit)
}

# S3 method for speedglm binomial family fit, takes BinDat data object:
logisfit.speedglmS3 <- function(BinDatObject, ...) {
  if (gvars$verbose) print("calling speedglm.wfit...")
  Xmat <- BinDatObject$getXmat
  Y_vals <- BinDatObject$getY
  # Xmat has 0 rows: return NA`s and avoid throwing exception
  if (nrow(Xmat) == 0L) {
    m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    m.fit <- try(speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial(), trace = FALSE), silent = TRUE)
    if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
      message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", m.fit)
      return(logisfit.glmS3(BinDatObject))
    }
  }
  fit <- list(coef = m.fit$coef, linkfun = "logit_linkinv", fitfunname = "speedglm", nobs = nrow(Xmat))
  if (gvars$verbose) print(fit$coef)
  class(fit) <- c(class(fit), c("speedglmS3"))
  return(fit)
}

# S3 method for h2o binomial family fit, takes BinDat data object:
logisfit.h2oglmS3 <- function(BinDatObject, ...) {
  if (gvars$verbose) print("calling h2o.glm...")
  # Xmat <- BinDatObject$getXmat
  # Y_vals <- BinDatObject$getY
  yname <- BinDatObject$outvar
  xnames <- BinDatObject$predvars
  subset_idx <- which(BinDatObject$subset_idx)
  if (length(subset_idx) == 0L) { # Xmat has 0 rows: return NA`s and avoid throwing exception
    m.fit <- list(coef = rep.int(NA_real_, length(xnames)))
  } else if (length(xnames) == 0L) {
    return(logisfit.speedglmS3(BinDatObject))
  } else {
    # Random Forests:
    # my.rf = h2o::h2o.randomForest(x = xnames, y = yname, training_frame = newH2Oframe, ntree = 100)

    # GBM:
    # for GBM to run need to make outcome into a factor:
    # newH2Oframe <- BinDatObject$DataStorageObject$H2O.dat.sVar[subset_idx,]
    # newH2Oframe[,yname] <- h2o::as.factor(newH2Oframe[,yname])
    # my.gbm <- h2o::h2o.gbm(x = xnames, y = yname, training_frame = newH2Oframe, distribution = "bernoulli")

    # GLM:
    m.fit <- try(h2o::h2o.glm(y = yname,
                              x = xnames,
                              intercept = TRUE,
                              training_frame = BinDatObject$DataStorageObject$H2O.dat.sVar[subset_idx,],
                              # training_frame = newH2Oframe,
                              family = "binomial",
                              standardize = TRUE,
                              solver = c("L_BFGS"), # solver = c("IRLSM"),
                              # remove_collinear_columns = TRUE,
                              max_iterations = 50,
                              lambda = 0L),
              silent = TRUE)

    # browser()
    # length(subset_idx)
    # print("h2o.glm fit"); print(m.fit)
    # print("h2o.glm coefficients"); print(m.fit@model$coefficients)
    # print("h2o.glm coefficients");
    # print(m.fit@parameters)
    # print(m.fit@allparameters)
    # print(m.fit@model)
    # str(m.fit@model)

    if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
      message("h2o::h2o.glm failed, falling back on speedglm; ", m.fit)
      return(logisfit.speedglmS3(BinDatObject))
    }
  }

  # assign the fitted coefficients in correct order (same as predictor order in xnames)
  out_coef <- vector(mode = "numeric", length = length(xnames)+1)
  out_coef[] <- NA
  names(out_coef) <- c("Intercept", xnames)
  out_coef[names(m.fit@model$coefficients)] <- m.fit@model$coefficients

  fit <- list(coef = out_coef, linkfun = "logit_linkinv", fitfunname = "h2o.glm", nobs = length(subset_idx), H2O.model.object = m.fit)
  if (gvars$verbose) print(fit$coef)
  class(fit) <- c(class(fit), c("h2oglmS3"))
  return(fit)
}

# Generic prediction fun for logistic regression coefs, predicts P(A = 1 | newXmat)
# No need for S3 for now, until need different pred. funs for different classes
# Does not handle cases with deterministic Anodes in the original data..
logispredict <- function(m.fit, BinDatObject) {
  assert_that(!is.null(BinDatObject$getXmat)); assert_that(!is.null(BinDatObject$subset_idx))
  # Set to default missing value for A[i] degenerate/degerministic/misval:
  # Alternative, set to default replacement val: pAout <- rep.int(gvars$misXreplace, newBinDatObject$n)
  pAout <- rep.int(gvars$misval, BinDatObject$n)
  if (sum(BinDatObject$subset_idx > 0)) {
    eta <- BinDatObject$getXmat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
    pAout[BinDatObject$subset_idx] <- match.fun(FUN = m.fit$linkfun)(eta)
  }
  return(pAout)
}

logispredict.long <- function(m.fit, BinDatObject) {
  assert_that(!is.null(BinDatObject$getXmat)); assert_that(!is.null(BinDatObject$subset_idx))
  assert_that(nrow(BinDatObject$getXmat)==length(private$Y_vals))
  pAout <- rep.int(gvars$misval, BinDatObject$n)
  if (sum(BinDatObject$subset_idx > 0)) {
    # -----------------------------------------------------------------
    # OBTAINING PREDICTIONS FOR LONG FORMAT P(Ind_j = 1 | Bin_j, W) BASED ON EXISTING POOLED FIT:
    # -----------------------------------------------------------------
    eta <- BinDatObject$getXmat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
    probA1 <- match.fun(FUN = m.fit$linkfun)(eta)
    # -----------------------------------------------------------------
    # GETTING ID-BASED PREDICTIONS (n) as cumprod of P(Ind_j = 1 | Bin_j, W) for j = 1, ..., K
    # -----------------------------------------------------------------
    ProbAeqa_long <- as.vector(probA1^(private$Y_vals) * (1L - probA1)^(1L - private$Y_vals))
    res_DT <- data.table(ID = BinDatObject$ID, ProbAeqa_long = ProbAeqa_long)
    res_DT <- res_DT[, list(cumprob = cumprod(ProbAeqa_long)), by = ID]
    data.table::setkeyv(res_DT, c("ID")) # sort by ID
    res_DT_short <- res_DT[unique(res_DT[, key(res_DT), with = FALSE]), mult = 'last']
    ProbAeqa <- res_DT_short[["cumprob"]]
    pAout[BinDatObject$subset_idx] <- ProbAeqa
  }
  return(pAout)
}


# Convert existing Bin matrix (Bin indicators) for continuous self$outvar into long format data.table with 3 columns:
# ID - row number; sVar_allB.j - bin indicators collapsed into one col; bin_ID - bin number identify for prev. columns
# automatically removed all missing (degenerate) bin indicators
binirized.to.DTlong <- function(BinsDat_wide, binID_seq, ID, bin_names, pooled_bin_name, name.sVar) {
  # Convert Bin matrix into a data.table (without data.frame as intermediate), with new row ID column:
  DT_BinsDat_wide <- data.table::as.data.table(BinsDat_wide)[, c("ID") := ID, with = FALSE]
  data.table::setcolorder(DT_BinsDat_wide, c("ID", names(DT_BinsDat_wide)[-ncol(DT_BinsDat_wide)]))
  # melt into long format:
  sVar_melt_DT <- melt(DT_BinsDat_wide,
                      id.vars = "ID",
                      measure.vars = bin_names,
                      value.name = pooled_bin_name,
                      variable.name = name.sVar,
                      variable.factor = FALSE,
                      na.rm = FALSE)
  nbin_rep <- rep(binID_seq, each = nrow(BinsDat_wide))
  # 1) Add bin_ID; 2) remove a column with Bin names; 3) remove all rows with NA value for outcome (degenerate bins)
  if (!is.data.table(sVar_melt_DT)) {
    class(sVar_melt_DT)
    stop("sVar_melt_DT is not a data.table")
  }
  sVar_melt_DT <- sVar_melt_DT[, c("bin_ID") := list(nbin_rep)][, name.sVar := NULL, with = FALSE][!is.na(get(pooled_bin_name))]
  data.table::setkeyv(sVar_melt_DT, c("ID", "bin_ID"))  # sort by ID, bin_ID to prepare for merge with predictors (sW)
  return(sVar_melt_DT)
}

# Prepare predictors (sW/X_mat) as data.table, adding row IDs for a join
# Join with sVar_melt_DT that is already in long format
# Need to check that created IDs match exactly for both datasets
join.Xmat = function(X_mat, sVar_melt_DT, ID) {
  nIDs <- length(unique(sVar_melt_DT[["ID"]]))
  assert_that(nIDs == nrow(X_mat))
  X_mat_DT <- data.table::as.data.table(X_mat)[, c("ID") := ID, with = FALSE]
  data.table::setkeyv(X_mat_DT, c("ID")) # sort by ID
  sVar_melt_DT <- sVar_melt_DT[X_mat_DT] # Merge long format (self$pooled_bin_name, binIDs) with predictors (sW)
  return(sVar_melt_DT)
}

## ---------------------------------------------------------------------
#' R6 class for storing the design matrix and the binary outcome for a single GLM (logistic) regression
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
#' \item{bin_names} - Names of the bins.
#' \item{ID} - Vector of observation IDs, \code{1:n}, used for pooling.
#' \item{pooled_bin_name} - Original name of the continuous covariate that was discretized into bins and then pooled.
#' \item{nbins} - Number of bins.
#' \item{outvar} - Outcome name.
#' \item{predvars} - Predictor names.
#' \item{pool_cont} - Perform pooling of bins?
#' \item{outvars_to_pool} - Outcome bin indicators to pool?
#' \item{subset_vars} - Defines the subset which would be used for fitting this model (logical, expression or indices).
#' \item{subset_idx} - Subset \code{subset_vars} converted to logical vector.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg)}}{Uses \code{reg} R6 \code{\link{RegressionClass}} object to instantiate a new storage container for a
#'   design matrix and binary outcome.}
#'   \item{\code{show()}}{ Print information on outcome and predictor names used in this regression model}
#'   \item{\code{newdata()}}{...}
#'   \item{\code{define.subset_idx(...)}}{...}
#'   \item{\code{setdata()}}{...}
#'   \item{\code{setdata.long()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{emptydata}}{...}
#'   \item{\code{emptyY}}{...}
#'   \item{\code{emptySubset_idx}}{...}
#'   \item{\code{emptyN}}{...}
#'   \item{\code{getXmat}}{...}
#'   \item{\code{getY}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.string is.flag
#' @export
BinDat <- R6Class(classname = "BinDat",
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    # reg = NULL,
    bin_names = NULL,
    ID = NULL,
    pooled_bin_name = NULL,
    # binID_seq = NULL,
    DataStorageObject = NULL,
    nbins = integer(),
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    pool_cont = logical(),
    outvars_to_pool = character(),
    ReplMisVal0 = logical(),
    n = NA_integer_,        # number of rows in the input data
    subset_vars = NULL,     # THE VAR NAMES WHICH WILL BE TESTED FOR MISSINGNESS AND WILL DEFINE SUBSETTING
    subset_expr = NULL,     # THE LOGICAL EXPRESSION (ONE) TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)

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

    setdata.long = function(data, ...) {
      assert_that(is.DataStorageClass(data))
      self$n <- data$nobs
      self$subset_idx <- self$define.subset_idx(data)
      if (!(data$active.bin.sVar %in% self$outvar)) { stop("currently binirized sVar does not match self$outvar argument") }

      # Setting up object fields related to pooling of continuous sA:
      self$pooled_bin_name <- data$pooled.bin.nm.sVar(self$outvar)
      self$bin_names <- self$outvars_to_pool

      if (gvars$verbose) {
        print("self$bin_names: "); print(self$bin_names)
        print("self$pooled_bin_name: "); print(self$pooled_bin_name)
        print("self$data$active.bin.sVar: "); print(self$data$active.bin.sVar)
        print("self$outvar: "); print(self$outvar)
        print("self$nbins: "); print(self$nbins)
      }

      binID_seq <- 1L:self$nbins
      BinsDat_wide <- data$get.dat.sVar(self$subset_idx, self$outvars_to_pool)
      self$ID <- as.integer(1:nrow(BinsDat_wide))

      # To grab bin Ind mat directly (prob a bit faster): BinsDat_wide <- data$dat.bin.sVar[self$subset_idx, ]
      BinsDat_long <- binirized.to.DTlong(BinsDat_wide = BinsDat_wide, binID_seq = binID_seq, ID = self$ID,
                                          bin_names = self$bin_names, pooled_bin_name = self$pooled_bin_name,
                                          name.sVar = self$outvar)
      sVar_melt_DT <- join.Xmat(X_mat = data$get.dat.sVar(self$subset_idx, self$predvars),
                                sVar_melt_DT = BinsDat_long, ID = self$ID)
      # prepare design matrix for modeling w/ glm.fit or speedglm.wfit:
      X_mat <- sVar_melt_DT[,c("bin_ID", self$predvars), with=FALSE][, c("Intercept") := 1] # select bin_ID + predictors, add intercept column
      setcolorder(X_mat, c("Intercept", "bin_ID", self$predvars)) # re-order columns by reference (no copy)
      self$ID <- sVar_melt_DT[["ID"]]
      private$X_mat <- as.matrix(X_mat)
      private$Y_vals <- sVar_melt_DT[, self$pooled_bin_name, with = FALSE][[1]] # outcome vector:

      if (gvars$verbose) {
        print("private$X_mat[1:10,]"); print(private$X_mat[1:10,])
        print("head(private$Y_vals)"); print(head(private$Y_vals, 100))
      }

      # **************************************
      # TO FINISH...
      # **************************************
      # if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
      #   private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
      #   colnames(private$X_mat) <- c("Intercept", self$predvars)
      # } else {
      #   # *** THIS IS THE ONLY LOCATION IN THE PACKAGE WHERE CALL TO DataStorageClass$get.dat.sVar() IS MADE ***
      #   private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.dat.sVar(self$subset_idx, self$predvars)))
      # To find and replace misvals in X_mat:
        if (self$ReplMisVal0) private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace
      # }
    }
  ),

  active = list( # 2 types of active bindings (w and wout args)
    emptydata = function() { private$X_mat <- NULL; self$DataStorageObject <- NULL },
    emptyY = function() { private$Y_vals <- NULL},
    emptySubset_idx = function() { self$subset_idx <- NULL },
    emptyN = function() { self$n <- NA_integer_ },
    getXmat = function() {private$X_mat},
    getY = function() {private$Y_vals}
  ),

  private = list(
    X_mat = NULL,
    Y_vals = NULL
  )
)