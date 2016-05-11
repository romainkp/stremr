#----------------------------------------------------------------------------------
# Classes for modelling regression models with binary outcome Bin ~ Xmat
#----------------------------------------------------------------------------------

logisfit <- function(datsum_obj) UseMethod("logisfit") # Generic for fitting the logistic model

# S3 method for glm binomial family fit, takes BinDat data object:
logisfit.glmS3 <- function(datsum_obj) {
  if (gvars$verbose) print("calling glm.fit...")
  Xmat <- datsum_obj$getXmat
  Y_vals <- datsum_obj$getY
    # Xmat has 0 rows: return NA's and avoid throwing exception:
  if (nrow(Xmat) == 0L) {
    m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    ctrl <- glm.control(trace = FALSE)
    # ctrl <- glm.control(trace = FALSE, maxit = 500)
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
logisfit.speedglmS3 <- function(datsum_obj) {
  if (gvars$verbose) print("calling speedglm.wfit...")
  Xmat <- datsum_obj$getXmat
  Y_vals <- datsum_obj$getY

  if (nrow(Xmat) == 0L) { # Xmat has 0 rows: return NA`s and avoid throwing exception
    m.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    # , maxit=1000
    m.fit <- try(speedglm::speedglm.wfit(X = Xmat, y = Y_vals, family = binomial(), trace = FALSE), silent = TRUE)
    if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
      message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", m.fit)
      return(logisfit.glmS3(datsum_obj))
    }
  }
  fit <- list(coef = m.fit$coef, linkfun = "logit_linkinv", fitfunname = "speedglm", nobs = nrow(Xmat))
  if (gvars$verbose) print(fit$coef)
  class(fit) <- c(class(fit), c("speedglmS3"))
  return(fit)
}

# S3 methods for getting coefs from fitted BinOutModel class object
coef.BinOutModel <- function(binoutmodel) {
  assert_that(binoutmodel$is.fitted)
  fit <- binoutmodel$getfit
  fit$coef
}

summary.BinOutModel <- function(binoutmodel) {
  assert_that(binoutmodel$is.fitted)
  fit <- binoutmodel$getfit
  append(list(reg = binoutmodel$show()), fit)
}

#' @import data.table
NULL

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
#' R6 class for storing the design matrix and binary outcome for a single logistic regression
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
#'   \item{\code{logispredict()}}{...}
#'   \item{\code{setdata.long()}}{...}
#'   \item{\code{logispredict.long()}}{...}
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
  # cloneable = FALSE,
  portable = TRUE,
  class = TRUE,
  public = list(
    # reg = NULL,
    bin_names = NULL,
    ID = NULL,
    pooled_bin_name = NULL,
    # binID_seq = NULL,
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
        # n_per_row <- 5
        # vec <- self$predvars
        # vec <- c(paste0(vec[-length(vec)], ", "), vec[length(vec)])
        # new_vec <- c(vec, rep("", length(vec) %/% n_per_row))
        # vec_byrow <- matrix(new_vec, ncol = n_per_row, byrow = TRUE)
        # new_pred_form <- paste0(as.vector(t(cbind(vec_byrow, c(rep(" \\", nrow(vec_byrow)-1), "")))), collapse="")
        # "P(" %+% self$outvar %+% "|" %+% paste(new_pred_form, collapse="") %+% ")" %+% "; \\
        #  Stratify: " %+% self$subset_expr
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
      # ******************************************************
      # NOTE: Below subsetting by call/expression is currently disabled, for speed & memory efficiency
      # all subsetting is done by subsetvars (variable name(s) must be non-missing)
      # ******************************************************
      } else if (is.call(self$subset_vars)) {
        # subset_idx <- data$evalsubst(subsetexpr = self$subset_vars)
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
    },

    logispredict.long = function(m.fit) {
      assert_that(!is.null(private$X_mat)); assert_that(!is.null(self$subset_idx))
      assert_that(nrow(private$X_mat)==length(private$Y_vals))
      pAout <- rep.int(gvars$misval, self$n)
      if (sum(self$subset_idx > 0)) {
        # -----------------------------------------------------------------
        # OBTAINING PREDICTIONS FOR LONG FORMAT P(Ind_j = 1 | Bin_j, W) BASED ON EXISTING POOLED FIT:
        # -----------------------------------------------------------------
        eta <- private$X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
        probA1 <- match.fun(FUN = m.fit$linkfun)(eta)
        # -----------------------------------------------------------------
        # GETTING ID-BASED PREDICTIONS (n) as cumprod of P(Ind_j = 1 | Bin_j, W) for j = 1, ..., K
        # -----------------------------------------------------------------
        ProbAeqa_long <- as.vector(probA1^(private$Y_vals) * (1L - probA1)^(1L - private$Y_vals))
        res_DT <- data.table(ID = self$ID, ProbAeqa_long = ProbAeqa_long)
        res_DT <- res_DT[, list(cumprob = cumprod(ProbAeqa_long)), by = ID]
        data.table::setkeyv(res_DT, c("ID")) # sort by ID
        res_DT_short <- res_DT[unique(res_DT[, key(res_DT), with = FALSE]), mult = 'last']
        ProbAeqa <- res_DT_short[["cumprob"]]
        # print("res_DT: "); print(res_DT)
        # print("res_DT w/ cumprob: "); print(res_DT)
        # print("res_DT_short: "); print(res_DT_short)
        # print("length(ProbAeqa): " %+% length(ProbAeqa))
        # print("head(ProbAeqa, 50)"); print(head(ProbAeqa, 50))
        pAout[self$subset_idx] <- ProbAeqa
      }
      return(pAout)
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

  active = list( # 2 types of active bindings (w and wout args)
    emptydata = function() { private$X_mat <- NULL },
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


## ---------------------------------------------------------------------
#' R6 class for fitting and making predictions for a single logistic regression with binary outcome B, P(B | PredVars)
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
#' \item{cont.sVar.flag} - Is the original outcome variable continuous?
#' \item{bw.j} - Bin width (interval length) for an outcome that is a bin indicator of a discretized continous outcome.
#' \item{glmfitclass} - Controls which package will be used for performing model fits (\code{glm} or \code{speedglm}).
#' \item{bindat} - Pointer to an instance of \code{BinDat} class that contains the data.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg)}}{Uses \code{reg} R6 \code{\link{RegressionClass}} object to instantiate a new model for a
#'   logistic regression with binary outcome.}
#'   \item{\code{show()}}{Print information on outcome and predictor names used in this regression model}
#'   \item{\code{fit()}}{...}
#'   \item{\code{copy.fit()}}{...}
#'   \item{\code{predict()}}{...}
#'   \item{\code{copy.predict()}}{...}
#'   \item{\code{predictAeqa()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{getoutvarnm}}{...}
#'   \item{\code{getoutvarval}}{...}
#'   \item{\code{getsubset}}{...}
#'   \item{\code{getprobA1}}{...}
#'   \item{\code{getfit}}{...}
#'   \item{\code{wipe.alldat}}{...}
#' }
#' @importFrom assertthat assert_that is.flag
#' @export
BinOutModel  <- R6Class(classname = "BinOutModel",
  # cloneable = FALSE,
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    outvar = character(),   # outcome name(s)
    predvars = character(), # names of predictor vars
    cont.sVar.flag = logical(),
    bw.j = numeric(),
    glmfitclass = "glmS3", # default glm fit class
    is.fitted = FALSE,
    bindat = NULL, # object of class BinDat that is used in fitting / prediction, never saved (need to be initialized with $new())

    initialize = function(reg, ...) {
      assert_that(is.flag(reg$useglm))
      if (!reg$useglm) self$glmfitclass <- "speedglmS3"

      self$outvar <- reg$outvar
      self$predvars <- reg$predvars

      self$bindat <- BinDat$new(reg = reg, ...) # postponed adding data in BinDat until self$fit() is called
      class(self$bindat) <- c(class(self$bindat), self$glmfitclass)
      if (gvars$verbose) {
        print("New instance of BinOutModel:"); print(self$show())
      }
      # Get the bin width (interval length) for the current bin name self$getoutvarnm (for discretized continuous sA only):
      self$cont.sVar.flag <- self$getoutvarnm %in% names(reg$intrvls.width)
      if (self$cont.sVar.flag) {
        intrvl.idx <- which(names(reg$intrvls.width) %in% self$getoutvarnm)
        if (length(intrvl.idx) > 1) stop("non-unique names for intrvls.width in RegressionClass")
        self$bw.j <- reg$intrvls.width[intrvl.idx]
      } else {
        self$bw.j <- 1L
      }
      invisible(self)
    },

    fit = function(overwrite = FALSE, data, ...) { # Move overwrite to a field? ... self$overwrite
      if (gvars$verbose) print("fitting the model: " %+% self$show())

      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of prev. fitted model unless explicitely asked
      self$bindat$newdata(newdata = data, ...) # populate bindat with X_mat & Y_vals
      private$m.fit <- logisfit(datsum_obj = self$bindat) # private$m.fit <- data_obj$logisfit or private$m.fit <- data_obj$logisfit()
      # alternative 2 is to apply data_obj method / method that fits the model
      self$is.fitted <- TRUE
      # **********************************************************************
      # to save RAM space when doing many stacked regressions no longer predicting in fit:
      # **********************************************************************
      # if (self$reg$pool_cont && length(self$reg$outvars_to_pool) > 1) {
      #   private$probAeqa <- self$bindat$logispredict.long(m.fit = private$m.fit)
      # } else {
      #   private$probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
      # }
      self$wipe.alldat
      invisible(self)
    },

    # take fitted BinOutModel class object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("BinOutModel" %in% class(bin.out.model))
      private$m.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # Predict the response P(Bin = 1|sW = sw);
    # uses private$m.fit to generate predictions for newdata:
    predict = function(newdata, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata)) {
        stop("must provide newdata for BinOutModel$predict()")
      }
      # re-populate bindat with new X_mat:
      self$bindat$newdata(newdata = newdata, getoutvar = FALSE, ...)
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        stop("BinOutModel$predict is not applicable to pooled regression, call BinOutModel$predictAeqa instead")
      } else {
        private$probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
      }
      self$bindat$emptydata  # Xmat in bindat is no longer needed, but subset, outvar & probA1 may be needed for private$probA1
      invisible(self)
    },

    # take BinOutModel class object that contains the predictions for P(A=1|sW) and save these predictions to self$
    copy.predict = function(bin.out.model) {
      assert_that("BinOutModel" %in% class(bin.out.model))
      assert_that(self$is.fitted)
      private$probA1 <- bin.out.model$getprobA1
    },

    # Predict the response P(Bin = b|sW = sw), which is returned invisibly;
    # Needs to know the values of b for prediction
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, bw.j.sA_diff) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      assert_that(!missing(newdata))
      self$bindat$newdata(newdata = newdata, getoutvar = TRUE) # populate bindat with new design matrix covars X_mat
      assert_that(is.logical(self$getsubset))
      n <- newdata$nobs
      # obtain predictions (likelihood) for response on fitted data (from long pooled regression):
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        probAeqa <- self$bindat$logispredict.long(m.fit = private$m.fit) # overwrite probA1 with new predictions:
      } else {
        # get predictions for P(sA[j]=1|sW=newdata) from newdata:
        probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
        indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
        assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
        probAeqa <- rep.int(1L, n) # for missing, the likelihood is always set to P(A = a) = 1.
        assert_that(!any(is.na(probA1[self$getsubset]))) # check that predictions P(A=1 | dmat) exist for all obs.
        probA1 <- probA1[self$getsubset]
        # discrete version for the joint density:
        probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
        # continuous version for the joint density:
        # probAeqa[self$getsubset] <- (probA1^indA) * exp(-probA1)^(1 - indA)
        # Alternative intergrating the last hazard chunk up to x:
        # difference of sA value and its left most bin cutoff: x - b_{j-1}
        if (!missing(bw.j.sA_diff)) {
          # + integrating the constant hazard all the way up to value of each sa:
          # probAeqa[self$getsubset] <- probAeqa[self$getsubset] * (1 - bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
          # cont. version of above:
          probAeqa[self$getsubset] <- probAeqa[self$getsubset] * exp(-bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
        }
      }
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      # private$probAeqa <- probAeqa # NOTE disabling internal saving of probAeqa
      # **********************************************************************
      return(probAeqa)
    },

    sampleA = function(newdata, bw.j.sA_diff) { # P(A^s[i]=a^s|W^s=w^s) - calculating the likelihood for indA[i] (n vector of a`s)
      assert_that(self$is.fitted)
      assert_that(!missing(newdata))
      # Don't want to subset by the outvar, since binarized mat for cat outcome is not re-created when just sampling
      # But need to reset it back when done
      temp_subset_vars <- self$bindat$subset_vars
      self$bindat$subset_vars <- self$bindat$subset_vars[!self$bindat$subset_vars %in% self$bindat$outvar]
      self$bindat$newdata(newdata = newdata, getoutvar = FALSE) # populate bindat with new design matrix covars X_mat

      assert_that(is.logical(self$getsubset))
      n <- newdata$nobs
      # obtain predictions (likelihood) for response on fitted data (from long pooled regression):
      if (self$bindat$pool_cont && length(self$bindat$outvars_to_pool) > 1) {
        stop("not implemented")
        # probAeqa <- self$bindat$logispredict.long(m.fit = private$m.fit) # overwrite probA1 with new predictions:
      } else {
        # get probability P(sA[j]=1|sW=newdata) from newdata, then sample from rbinom
        probA1 <- self$bindat$logispredict(m.fit = private$m.fit)
        sampleA <- rep.int(0L, n)
        sampleA[self$getsubset] <- rbinom(n = n, size = 1, prob = probA1)

      #   indA <- newdata$get.outvar(self$getsubset, self$getoutvarnm) # Always a vector of 0/1
      #   assert_that(is.integerish(indA)) # check that obsdat.sA is always a vector of of integers
      #   probAeqa <- rep.int(1L, n) # for missing, the likelihood is always set to P(A = a) = 1.
      #   assert_that(!any(is.na(probA1[self$getsubset]))) # check that predictions P(A=1 | dmat) exist for all obs.
      #   probA1 <- probA1[self$getsubset]
      #   # discrete version for the joint density:
      #   probAeqa[self$getsubset] <- probA1^(indA) * (1 - probA1)^(1L - indA)
      #   # continuous version for the joint density:
      #   # probAeqa[self$getsubset] <- (probA1^indA) * exp(-probA1)^(1 - indA)
      #   # Alternative intergrating the last hazard chunk up to x:
      #   # difference of sA value and its left most bin cutoff: x - b_{j-1}
      #   if (!missing(bw.j.sA_diff)) {
      #     # + integrating the constant hazard all the way up to value of each sa:
      #     # probAeqa[self$getsubset] <- probAeqa[self$getsubset] * (1 - bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
      #     # cont. version of above:
      #     probAeqa[self$getsubset] <- probAeqa[self$getsubset] * exp(-bw.j.sA_diff[self$getsubset]*(1/self$bw.j)*probA1)^(indA)
      #   }
      }
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      self$wipe.alldat
      self$bindat$subset_vars <- temp_subset_vars
      # private$probAeqa <- probAeqa # NOTE disabling internal saving of probAeqa
      # **********************************************************************
      return(sampleA)
    },

    get.fits = function(format_table = FALSE) {
      coef_out <- private$m.fit$coef
      nobs <- private$m.fit$nobs
      if (format_table) {
        coef_out <- data.frame(Coef = coef_out)
        coef_out <- cbind(names(private$m.fit$coef), coef_out)
        rownames(coef_out) <- NULL
        colnames(coef_out) <- c("Terms", "Coefficients")
      }
      return(list(c(self$show(print_format = FALSE), list(nobs = nobs, coef = coef_out))))
    },
    show = function(print_format = TRUE) {self$bindat$show(print_format)}
    # ,
    # # return new R6 object that only contains a copy of the fits in self
    # clone = function(deep = TRUE) {
    #   BinOutModel$new(reg = reg, ...)
    #   assert_that("BinOutModel" %in% class(bin.out.model))
    #   assert_that(self$is.fitted)
    #   private$probA1 <- bin.out.model$getprobA1
    # }
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