#-----------------------------------------------------------------------------
# DataStorageClass CLASS STRUCTURE:
#-----------------------------------------------------------------------------
# Creates and stores the combined summary measure matrix (sW,sA) in DataStorageClass class;
# Contains Methods for:
  # *) detecting sVar types (detect.col.types);
  # *) defining interval cuttoffs for continuous sVar (define.intervals)
  # *) turning continuous sVar into categorical (discretize.sVar)
  # *) creating binary indicator matrix for continous/categorical sVar (binirize.sVar, binirize.cat.sVar)
  # *) creating design matrix (Xmat) based on predvars and row subsets (evalsubst)

## ---------------------------------------------------------------------
# Detecting vector types: sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
## ---------------------------------------------------------------------
detect.col.types <- function(sVar_mat){
  detect_vec_type <- function(vec) {
    vec_nomiss <- vec[!gvars$misfun(vec)]
    nvals <- data.table::uniqueN(vec_nomiss)
    if (nvals <= 2L) {
      sVartypes$bin
    } else if ((nvals <= maxncats) && (is.integerish(vec_nomiss))) {
      sVartypes$cat
    } else {
      sVartypes$cont
    }
  }
  assert_that(is.integerish(getopt("maxncats")) && getopt("maxncats") > 1)
  maxncats <- getopt("maxncats")
  sVartypes <- gvars$sVartypes
  if (is.matrix(sVar_mat)) { # for matrix:
    return(as.list(apply(sVar_mat, 2, detect_vec_type)))
  } else if (is.data.table(sVar_mat)) { # for data.table:
    return(as.list(sVar_mat[, lapply(.SD, detect_vec_type)]))
  } else {
    stop("unrecognized sVar_mat class: " %+% class(sVar_mat))
  }
}

## ---------------------------------------------------------------------
# Normalizing / Defining bin intervals / Converting contin. to ordinal / Converting ordinal to bin indicators
## ---------------------------------------------------------------------
normalize <- function(x) {
  if (abs(max(x) - min(x)) > gvars$tolerr) { # Normalize to 0-1 only when x is not constant
    return((x - min(x)) / (max(x) - min(x)))
  } else {  # What is the thing to do when x constant? Set to abs(x), abs(x)/x or 0???
    return(x)
  }
}
# Define bin cutt-offs for continuous x:
define.intervals <- function(x, nbins, bin_bymass, bin_bydhist, max_nperbin) {
  x <- x[!gvars$misfun(x)]  # remove missing vals
  nvals <- length(unique(x))
  if (is.na(nbins)) nbins <- as.integer(length(x) / max_nperbin)
  # if nbins is too high, for ordinal, set nbins to n unique obs and cancel quantile based interval defns
  if (nvals < nbins) {
    nbins <- nvals
    bin_bymass <- FALSE
  }
  if (abs(max(x) - min(x)) > gvars$tolerr) {  # when x is not constant
    if ((bin_bymass) & !is.null(max_nperbin)) {
      if ((length(x) / max_nperbin) > nbins) nbins <- as.integer(length(x) / max_nperbin)
    }
    intvec <- seq.int(from = min(x), to = max(x) + 1, length.out = (nbins + 1)) # interval type 1: bin x by equal length intervals of 0-1
  } else {  # when x is constant, force the smallest possible interval to be at least [0,1]
    intvec <- seq.int(from = min(0L, min(x)), to = max(1L, max(x)), length.out = (nbins + 1))
  }
  if (bin_bymass) {
    intvec <- quantile(x = x, probs = normalize(intvec)) # interval type 2: bin x by mass (quantiles of 0-1 intvec as probs)
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  } else if (bin_bydhist) {
    stop("... binning continuous variable: dhist bin definitions are no longer available...")
    # intvec <- dhist(x, plot = FALSE, nbins = nbins)$xbr
    # intvec[1] <- intvec[1] - 0.01
    # intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  }
  # adding -Inf & +Inf as leftmost & rightmost cutoff points to make sure all future data points end up in one of the intervals:
  intvec <- c(min(intvec)-1000L, intvec, max(intvec)+1000L)
  return(intvec)
}
# Turn any x into ordinal (1, 2, 3, ..., nbins) for a given interval cutoffs (length(intervals)=nbins+1)
make.ordinal <- function(x, intervals) findInterval(x = x, vec = intervals, rightmost.closed = TRUE)
# Make dummy indicators for ordinal x (sA[j]) Approach used: creates B_j that jumps to 1 only once and stays 1 (degenerate) excludes reference category (last)
make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms, levels = 1:nbins) {
  n <- length(x.ordinal)
  new.cats <- 1:nbins
  dummies_mat <- matrix(1L, nrow = n, ncol = length(new.cats))
  for(cat in new.cats[-length(new.cats)]) {
    subset_Bj0 <- x.ordinal > levels[cat]
    dummies_mat[subset_Bj0, cat] <- 0L
    subset_Bjmiss <- x.ordinal < levels[cat]
    dummies_mat[subset_Bjmiss, cat] <- gvars$misval
  }
  dummies_mat[, new.cats[length(new.cats)]] <- gvars$misval
  colnames(dummies_mat) <- bin.nms
  return(dummies_mat)
}

# put the first level (category) as a reference (last category) for categorical censoring
make.bins_mtx_cens <- function(x.ordinal, nbins, bin.nms, levels = 1:nbins, ref.level = gvars$noCENScat) {
  n <- length(x.ordinal)
  new.ref.level <- levels[length(levels)]+1
  levels <- c(levels[-which(levels %in% ref.level)], new.ref.level)
  x.ordinal[as.integer(x.ordinal)==as.integer(ref.level)] <- new.ref.level

  new.cats <- 1:nbins
  dummies_mat <- matrix(1L, nrow = n, ncol = length(new.cats))
  for(cat in new.cats[-length(new.cats)]) {
    subset_Bj0 <- x.ordinal > levels[cat]
    dummies_mat[subset_Bj0, cat] <- 0L
    subset_Bjmiss <- x.ordinal < levels[cat]
    dummies_mat[subset_Bjmiss, cat] <- gvars$misval
  }
  dummies_mat[, new.cats[length(new.cats)]] <- gvars$misval
  colnames(dummies_mat) <- bin.nms
  return(dummies_mat)
}

is.H2OFrame <- function(fr)  base::`&&`(!missing(fr), class(fr)[1]=="H2OFrame")

## ---------------------------------------------------------------------
#' R6 class for storing, managing, subsetting and manipulating the input data.
#'
#'  The class \code{DataStorageClass} is the only way the package uses to access the input data.
#'  The evaluated summary measures from sVar.object are stored as a matrix (\code{private$.mat.sVar}).
#'  Contains methods for replacing missing values with default in \code{gvars$misXreplace}.
#'  Also contains method for detecting / setting sVar variable type (binary, categor, contin).
#'  Contains methods for combining, subsetting, discretizing & binirizing summary measures \code{(sW,sA)}.
#'  For continous sVar this class provides methods for detecting / setting bin intervals,
#'  normalization, disretization and construction of bin indicators.
#'  The pointers to this class get passed on to \code{GenericModel} functions: \code{$fit()},
#'  \code{$predict()} and \code{$predictAeqa()}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#'    \item{\code{YnodeVals}}
#'    \item{\code{det.Y}}
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(Odata, nodes, YnodeVals, det.Y, ...)}}{...}
#'   \item{\code{def.types.sVar(type.sVar = NULL)}}{...}
#'   \item{\code{fixmiss_sVar()}}{...}
#'   \item{\code{set.sVar(name.sVar, new.type)}}{...}
#'   \item{\code{set.sVar.type(name.sVar, new.type)}}{...}
#'   \item{\code{get.sVar(name.sVar, new.sVarVal)}}{...}
#'   \item{\code{replaceOneAnode(AnodeName, newAnodeVal)}}{...}
#'   \item{\code{replaceManyAnodes(Anodes, newAnodesMat)}}{...}
#'   \item{\code{addYnode(YnodeVals, det.Y)}}{...}
#'   \item{\code{evalsubst(subset_exprs, subset_vars)}}{...}
#'   \item{\code{get.dat.sVar(rowsubset = TRUE, covars)}}{...}
#'   \item{\code{get.outvar(rowsubset = TRUE, var)}}{...}
#'   \item{\code{bin.nms.sVar(name.sVar, nbins)}}{...}
#   \item{\code{pooled.bin.nm.sVar(name.sVar)}}{...}
#'   \item{\code{detect.sVar.intrvls(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin)}}{...}
#'   \item{\code{detect.cat.sVar.levels(name.sVar)}}{...}
#'   \item{\code{binirize.sVar(name.sVar, ...)}}{...}
#'   \item{\code{get.sVar.bw(name.sVar, intervals)}}{...}
#'   \item{\code{get.sVar.bwdiff(name.sVar, intervals)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{nobs}}{...}
#'    \item{\code{ncols.sVar}}{...}
#'    \item{\code{names.sVar}}{...}
#'    \item{\code{type.sVar}}{Named list of length \code{ncol(private$.mat.sVar)} with \code{sVar} variable types: "binary"/"categor"/"contin".}
#'    \item{\code{dat.sVar}}{...}
#'    \item{\code{ord.sVar}}{Ordinal (categorical) transformation of a continous covariate \code{sVar}.}
#'    \item{\code{active.bin.sVar}}{Name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in private$.mat.bin.sVar)}
#'    \item{\code{dat.bin.sVar}}{...}
#'    \item{\code{emptydat.sVar}}{...}
#'    \item{\code{emptydat.bin.sVar}}{...}
#'    \item{\code{noNA.Ynodevals}}{...}
#'    \item{\code{nodes}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.flag
#' @export
DataStorageClass <- R6Class(classname = "DataStorageClass",
  portable = TRUE,
  class = TRUE,
  public = list(
    Qlearn.fit = NULL,
    interventionNodes.g0 = NULL,
    interventionNodes.gstar = NULL,
    gstarNodes_stoch = NULL,
    modelfits.g0 = NULL,
    modelfit.gC = NULL,
    modelfit.gA = NULL,
    modelfit.gN = NULL,
    new.factor.names = NULL,
    noCENScat = 0L,         # The level (integer) that indicates CONTINUATION OF FOLLOW-UP for ALL censoring variables
    YnodeVals = NULL,       # Values of the binary outcome (Ynode) in observed data where det.Y = TRUE obs are set to NA
    det.Y = NULL,           # Logical vector, where YnodeVals[det.Y==TRUE] are deterministic (0 or 1)
    curr_data_A_g0 = TRUE,  # is the current data in OdataDT generated under observed (g0)? If FALSE, current data is under g.star (intervention)
    fold_column = NULL,
    nfolds = NULL,

    initialize = function(Odata, nodes, YnodeVals, det.Y, noCENScat,...) {
      assert_that(is.data.frame(Odata) | is.data.table(Odata))
      self$curr_data_A_g0 <- TRUE
      self$dat.sVar <- data.table(Odata) # makes a copy of the input data (shallow)
      # alternative is to set it without copying Odata
      # setDT(Odata); self$dat.sVar <- Odata

      # set the keys for quick search!
      setkeyv(self$dat.sVar, cols = c(nodes$IDnode, nodes$tnode))

      if (!missing(noCENScat)) {
        self$noCENScat <- noCENScat
      } else {
        self$noCENScat <- gvars$noCENScat
      }

      if (!missing(nodes)) self$nodes <- nodes

      if (!missing(YnodeVals)) self$addYnode(YnodeVals = YnodeVals, det.Y = det.Y)

      if (gvars$verbose) print("...detecting the type of each input column...")
      self$def.types.sVar() # Define the type of each sVar[i]: bin, cat or cont

      invisible(self)
    },

    # add protected Y nodes to private field and set to NA all determinisitc Y values for public field YnodeVals
    addYnode = function(YnodeVals, det.Y) {
        if (missing(det.Y)) det.Y <- rep.int(FALSE, length(YnodeVals))
        self$noNA.Ynodevals <- YnodeVals  # Adding actual observed Y as protected (without NAs)
        self$YnodeVals <- YnodeVals
        self$YnodeVals[det.Y] <- NA       # Adding public YnodeVals & setting det.Y values to NA
        self$det.Y <- det.Y
    },

    # ---------------------------------------------------------------------
    # Eval the subsetting expression (in the environment of the data.table "data" + global constants "gvars"):
    # ---------------------------------------------------------------------
    # Could also do evaluation in a special env with a custom subsetting fun '[' that will dynamically find the correct dataset that contains
    # sVar.name (dat.sVar or dat.bin.sVar) and will return sVar vector
    evalsubst = function(subset_vars, subset_exprs = NULL) {
      res <- rep.int(TRUE, self$nobs)
      if (!missing(subset_vars)) {
        assert_that(is.character(subset_vars))
        for (subsetvar in subset_vars) {
          # *) find the var of interest (in self$dat.sVar or self$dat.bin.sVar), give error if not found
          sVar.vec <- self$get.outvar(var = subsetvar)
          assert_that(!is.null(sVar.vec))
          # *) reconstruct correct expression that tests for missing values
          res <- res & (!gvars$misfun(sVar.vec))
        }
      }
      if (!is.null(subset_exprs)) {
        if (is.logical(subset_exprs)) {
          res <- res & subset_exprs
        } else if (is.character(subset_exprs)){
          # ******************************************************
          # data.table evaluation of the logical subset expression
          # Note: This can be made a lot more faster by also keying data.table on variables in eval(parse(text = subset_exprs))
          # ******************************************************
          res.tmp <- self$dat.sVar[, eval(parse(text = subset_exprs)), by = get(self$nodes$ID)][["V1"]]
          assert_that(is.logical(res.tmp))
          res <- res & res.tmp
          # ******************************************************
          # THIS WAS A BOTTLENECK: for 500K w/ 1000 bins: 4-5sec
          # REPLACING WITH env that is made of data.frames instead of matrices
          # ******************************************************
          # eval.env <- c(data.frame(self$dat.sVar), data.frame(self$dat.bin.sVar), as.list(gvars))
          # res <- try(eval(subset_exprs, envir = eval.env, enclos = baseenv())) # to evaluate vars not found in data in baseenv()
          # self$dat.sVar[eval(),]
          # stop("disabled for memory/speed efficiency")
        }
      }
      return(res)
    },

    # ---------------------------------------------------------------------
    # Functions for subsetting/returning covariate design mat for BinaryOutcomeModel Class or outcome variable
    # ---------------------------------------------------------------------
    get.dat.sVar = function(rowsubset = TRUE, covars) {
      if (!missing(covars)) {
        if (length(unique(colnames(self$dat.sVar))) < length(colnames(self$dat.sVar))) {
          warning("repeating column names in the final data set; please check for duplicate summary measure / node names")
        }
        # columns to select from main design matrix (in the same order as listed in covars):
        sel.sWsA <- intersect(covars, colnames(self$dat.sVar))
        if (is.matrix(self$dat.sVar)) {
          dfsel <- self$dat.sVar[rowsubset, sel.sWsA, drop = FALSE] # data stored as matrix
        } else if (is.data.table(self$dat.sVar)) {
          dfsel <- self$dat.sVar[rowsubset, sel.sWsA, drop = FALSE, with = FALSE] # data stored as data.table
        } else {
          stop("self$dat.sVar is of unrecognized class: " %+% class(self$dat.sVar))
        }
        # columns to select from binned continuous/cat var matrix (if it was previously constructed):
        if (!is.null(self$dat.bin.sVar)) {
          sel.binsA <- intersect(covars, colnames(self$dat.bin.sVar))
        } else {
          sel.binsA <- NULL
        }
        if (length(sel.binsA)>0) {
          dfsel <- cbind(dfsel, self$dat.bin.sVar[rowsubset, sel.binsA, drop = FALSE])
        }
        found_vars <- covars %in% colnames(dfsel)
        if (!all(found_vars)) stop("some covariates can't be found: "%+%
                                    paste(covars[!found_vars], collapse=","))
        return(dfsel)
      } else {
        return(self$dat.sVar[rowsubset, , drop = FALSE])
      }
    },

    get.outvar = function(rowsubset = TRUE, var) {
      if (length(self$nodes) < 1) stop("DataStorageClass$nodes list is empty!")
      if (var %in% self$names.sVar) {
        out <- self$dat.sVar[rowsubset, var, with = FALSE]
      } else if (var %in% colnames(self$dat.bin.sVar)) {
        out <- self$dat.bin.sVar[rowsubset, var]
      } else if ((var %in% self$nodes$Ynode) && !is.null(self$YnodeVals)) {
        out <- self$YnodeVals[rowsubset]
      } else {
        stop("requested variable " %+% var %+% " does not exist in DataStorageClass!")
      }
      if ((is.list(out) || is.data.table(out)) && (length(out)>1)) {
        stop("selecting regression outcome covariate resulted in more than one column: " %+% var)
      } else if (is.list(out) || is.data.table(out)) {
        return(out[[1]])
      } else {
        return(out)
      }
    },

    # --------------------------------------------------
    # Create a matrix of dummy bin indicators for categorical/continuous sVar
    # --------------------------------------------------
    binirize.sVar = function(name.sVar, ...) {
      private$.active.bin.sVar <- name.sVar
      if (self$is.sVar.cont(name.sVar)) {
        private$binirize.cont.sVar(name.sVar, ...)
      } else if (self$is.sVar.cat(name.sVar)) {
        if (self$is.sVar.CENS(name.sVar)) {
          private$binirize.cat.CENS.sVar(name.sVar, ...)
        } else {
          private$binirize.cat.sVar(name.sVar, ...)
        }
      } else {
        stop("...can only call $binirize.sVar for continuous or categorical sVars...")
      }
    },

    # ------------------------------------------------------------------------------------------------------------
    # Binning methods for categorical/continuous sVar
    # ------------------------------------------------------------------------------------------------------------
    bin.nms.sVar = function(name.sVar, nbins) { name.sVar%+%"_"%+%"B."%+%(1L:nbins) }, # Return names of bin indicators for sVar:
    # pooled.bin.nm.sVar = function(name.sVar) { name.sVar %+% "_allB.j" },
    detect.sVar.intrvls = function(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin) {
      tol.int <- 0.001
      int <- define.intervals(x = self$get.sVar(name.sVar), nbins = nbins, bin_bymass = bin_bymass, bin_bydhist = bin_bydhist, max_nperbin = max_nperbin)
      diffvec <- diff(int)
      if (sum(abs(diffvec) < tol.int) > 0) {
        if (gvars$verbose) {
          message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+%
                  (length(int)-1) %+% " to " %+% (length(int[diffvec >= tol.int])-1) %+% " due to too few obs.")
          print("old intervals: "); print(as.numeric(int))
        }
        # Just taking unique interval values is insufficient
        # Instead need to drop all intervals that are "too close" to each other based on some tol value
        # remove all intervals (a,b) where |b-a| < tol.int, but always keep the very first interval (int[1])
        int <- c(int[1], int[2:length(int)][abs(diffvec) >= tol.int])
        if (gvars$verbose) print("new intervals: "); print(as.numeric(int))
      }
      return(int)
    },
    detect.cat.sVar.levels = function(name.sVar) {
      levels <- sort(unique(self$get.sVar(name.sVar)))
      return(levels)
    },
    # return the bin widths vector for the discretized continuous sVar (private$.ord.sVar):
    get.sVar.bw = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      intrvls.width <- diff(intervals)
      intrvls.width[intrvls.width <= gvars$tolerr] <- 1
      ord.sVar_bw <- intrvls.width[self$ord.sVar]
      return(ord.sVar_bw)
    },
   # return the bin widths vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bwdiff = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      ord.sVar_leftint <- intervals[self$ord.sVar]
      diff_bw <- self$get.sVar(name.sVar) - ord.sVar_leftint
      return(diff_bw)
    },

    # --------------------------------------------------
    # Replace all missing (NA) values with a default integer (0)
    # --------------------------------------------------
    fixmiss_sVar = function() {
      if (is.matrix(self$dat.sVar)) {
        private$fixmiss_sVar_mat()
      } else if (is.data.table(self$dat.sVar)) {
        private$fixmiss_sVar_DT()
      } else {
        stop("self$dat.sVar is of unrecognized class")
      }
    },

    # --------------------------------------------------
    # Methods for sVar types. Define the type (class) of each variable (sVar) in input data: gvars$sVartypes$bin,  gvars$sVartypes$cat or gvars$sVartypes$cont
    # --------------------------------------------------
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar), otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    def.types.sVar = function(type.sVar = NULL) {
      if (is.null(type.sVar)) {
        private$.type.sVar <- detect.col.types(self$dat.sVar)
      } else {
        n.sVar <- length(self$names.sVar)
        len <- length(type.sVar)
        assert_that((len == n.sVar) || (len == 1L))
        if (len == n.sVar) { # set types for each variable
          assert_that(is.list(type.sVar))
          assert_that(all(names(type.sVar) %in% self$names.sVar))
        } else { # set one type for all vars
          assert_that(is.string(type.sVar))
          type.sVar <- as.list(rep(type.sVar, n.sVar))
          names(type.sVar) <- self$names.sVar
        }
        private$.type.sVar <- type.sVar
      }
      invisible(self)
    },
    set.sVar.type = function(name.sVar, new.type) { private$.type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { private$.type.sVar } else { private$.type.sVar[[name.sVar]] } },
    is.sVar.bin = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$bin },
    is.sVar.cat = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$cat },
    is.sVar.cont = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$cont },
    is.sVar.CENS = function(name.sVar) { name.sVar %in% self$nodes$Cnodes },

    # ---------------------------------------------------------------------
    # Directly replace variable(s) in the storage data.table (by reference)
    # ---------------------------------------------------------------------
    get.sVar = function(name.sVar) {
      x <- self$dat.sVar[, name.sVar, with=FALSE]
      if (is.list(x) || is.data.table(x) || is.data.frame(x)) x <- x[[1]]
      return(x)
    },
    set.sVar = function(name.sVar, new.sVarVal) {
      # assert_that(is.integer(new.sVarVal) | is.numeric(new.sVarVal))
      assert_that(length(new.sVarVal)==self$nobs | length(new.sVarVal)==1)
      assert_that(name.sVar %in% colnames(self$dat.sVar))
      self$dat.sVar[, (name.sVar) := new.sVarVal]
      invisible(self)
    },

    # create a back-up of the observed input gstar nodes (created by user in input data):
    # needs to know how to add new columns (not backed up yet) TO SAME backup data.table
    backupNodes = function(nodes) {
      nodes <- nodes[!is.null(nodes)]
      for (node in nodes) CheckVarNameExists(self$dat.sVar, node)
      private$.saveGstarsDT <- self$dat.sVar[, nodes, with = FALSE]
      invisible(return(self))
    },

    restoreNodes = function(nodes) {
      nodes <- nodes[!is.null(nodes)]
      for (node in nodes) CheckVarNameExists(self$dat.sVar, node)
      if (is.null(private$.saveGstarsDT)) stop("Nodes in dat.sVar cannot be restored, private$.saveGstarsDT is null!")
      self$dat.sVar[, (nodes) := private$.saveGstarsDT[, nodes, with = FALSE]]
      # self$dat.sVar[, (nodes) := private$.saveGstarsDT[, nodes, with = FALSE], with = FALSE]
      invisible(return(self))
    },

    # Modify the values in node nodes_to_repl in self$dat.sVar with values from source_for_repl using only the IDs in subset_idx:
    replaceNodesVals = function(subset_idx, nodes_to_repl = intervened_NODE, source_for_repl = NodeNames) {
      for (node_idx in seq_along(nodes_to_repl)) {
        if (length(subset_idx) > 0) {
        # if (sum(subset_idx) > 0) {
          source_node <- self$dat.sVar[subset_idx, (source_for_repl[node_idx]), with = FALSE][[source_for_repl[node_idx]]]
          self$dat.sVar[subset_idx, (nodes_to_repl[node_idx]) := as.numeric(source_node)]
        }
      }
      invisible(return(self))
    },

    # ## ---------------------------------------------------------------------
    # # Swap re-saved Anodes and summaries sA with those in main data.table
    # ## ---------------------------------------------------------------------
    # swapAnodes = function(Anodes) {
    #   if (missing(Anodes)) Anodes <- self$nodes$Anodes
    #   # 1) Save the current values of Anodes and sA in the data:
    #   temp.Anodes <- self$dat.sVar[, Anodes, with = FALSE]
    #   if (!is.null(private$.sA_g0_DT) && !is.null(self$save_sA_Vars)) {
    #     temp.sA <- self$dat.sVar[, self$save_sA_Vars, with = FALSE]
    #   } else {
    #     temp.sA <- NULL
    #   }
    #   # 2) Restore previously saved Anodes / sA into the data:
    #   self$restoreAnodes()
    #   # 3) Over-write the back-up values with new ones:
    #   private$.A_g0_DT <- temp.Anodes
    #   private$.sA_g0_DT <- temp.sA
    #   # 4) Reverse the indicator of current data Anodes:
    #   self$curr_data_A_g0 <- !self$curr_data_A_g0
    #   invisible(self)
    # },
    # replaceOneNode = function(NodeName, newNodeVal) {
    #   self$set.sVar(NodeName, newNodeVal)
    #   invisible(self)
    # },
    # replaceManyNodes = function(Nodes, newNodesMat) {
    #   assert_that(is.matrix(newNodesMat))
    #   assert_that(ncol(newNodesMat) == length(Nodes))
    #   for (col in Nodes) {
    #     idx <- which(Nodes %in% col)
    #     assert_that(col %in% colnames(newNodesMat))
    #     assert_that(col %in% colnames(self$dat.sVar))
    #     self$dat.sVar[, (col) := newNodesMat[, idx]]
    #   }
    #   invisible(self)
    # },
    # # Replace a column or columns in private$Xmat with new values
    # replaceCols = function(subset_idx, colnames) {
    #   for (colname in colnames) {
    #     private$Xmat[, colname]  <- self$get.dat.sVar(subset_idx, colname)
    #   }
    #   return(invisible(self))
    # },

    # swap node names in the data.table: current -> target and target -> current
    swapNodes = function(current, target) {
      # if current and target have the same node names, will result in error, so exclude
      common_names <- intersect(current, target)
      current <- current[!current %in% common_names]
      target <- target[!target %in% common_names]

      data.table::setnames(self$dat.sVar, old = current, new = current %+% ".temp.current")
      data.table::setnames(self$dat.sVar, old = target, new = current)
      data.table::setnames(self$dat.sVar, old = current %+% ".temp.current", new = target)
    },

    eval_rule_followers = function(NodeName, gstar.NodeName) {
      dat.rulefollow <- self$dat.sVar[, c(self$nodes$IDnode, NodeName, gstar.NodeName), with = FALSE]
      setkeyv(dat.rulefollow, cols = self$nodes$IDnode)
      # likelihood of N(t) under counterfactual g^*:
      dat.rulefollow[, "rule.follower.byt" := get(gstar.NodeName)^get(NodeName) * (1L-get(gstar.NodeName))^(1-get(NodeName))]
      # rule follower whenever the probability of observing the current history \bar{(O(t))} is > 0
      rule_followers_idx <- dat.rulefollow[, "rule.followers" := as.integer(cumprod(rule.follower.byt) > 0), by = eval(self$nodes$IDnode)][["rule.followers"]]
      return(as.logical(rule_followers_idx))
    },

    eval_uncensored = function() {
      return(self$dat.sVar[, list(uncensored_idx = as.logical(rowSums(.SD, na.rm = TRUE) == eval(self$noCENScat))), .SDcols = self$nodes$Cnodes][["uncensored_idx"]])
    },

    define.stoch.nodes = function(NodeNames) {
      gstarNodes_stoch <- vector(mode = "logical", length = length(NodeNames))
      names(gstarNodes_stoch) <- NodeNames
      for (gstarNode in NodeNames) {
        gstarNodes_stoch[[gstarNode]] <- !is.integerish(self$get.sVar(gstarNode)[!is.na(self$get.sVar(gstarNode))])
      }
      self$gstarNodes_stoch <- gstarNodes_stoch
      return(gstarNodes_stoch)
    },

    check_norows_after_event = function() {
      rows.exist.after.FAIL <- self$dat.sVar[, {outrow = which(get(self$nodes$Ynode) %in% 1L);
                                                LastRowIdx = .I[.N];
                                                StartRemoveRow = LastRowIdx - (.N - outrow) + 1;
                                                rows.exist =  (outrow != .N) && length(outrow) > 0;
                                                list(rows.exist = rows.exist, StartRemoveRow = StartRemoveRow, LastRowIdx = LastRowIdx)},
                                              by = eval(self$nodes$IDnode)][!is.na(rows.exist), ]

      if (any(rows.exist.after.FAIL[["rows.exist"]])) {
        rows.exist.after.FAIL <- rows.exist.after.FAIL[rows.exist.after.FAIL[["rows.exist"]], ]
        rows.exist.after.FAIL[, removeIDX := list(list((StartRemoveRow:LastRowIdx))), by = eval(self$nodes$IDnode)]
        idx.to.remove <- unlist(rows.exist.after.FAIL[["removeIDX"]])
        self$dat.sVar <- self$dat.sVar[!idx.to.remove, ]
        msg <- "Found " %+% length(idx.to.remove) %+% " extra rows after time-to-event (outcome) occurrence.
        Such rows aren't allowed and thus have been removed automatically.
        Unfortunately, this requires making a copy of the input data itself, doublying the amount of used RAM.
        If available RAM becomes an issue, consider first removing the original input data from memory by typing 'rm(...)',
        since its no longer needed for running stremr.
        To conserve RAM please consider removing these extra rows from your data!"
        warning(msg)
        message(msg)
      }
      invisible(return(self))
    },

    # ---------------------------------------------------------------------------
    # Cast long format data into wide format:
    # bslcovars - names of covariates that shouldn't be cast (remain invariant with t)
    # TO DO - add excludevars arg to exclude covariates
    # ---------------------------------------------------------------------------
    convert.to.wide = function(bslcovars) {
      # excludevars
      # dt = rbind(data.table(ID=1, x=sample(5,20,TRUE), y = sample(5,20,TRUE), t=1:20), data.table(ID=2, x=sample(5,15,TRUE), y = sample(5,15,TRUE), t=1:15))
      # dcast(dt, formula="ID ~ t", value.var=c("x", "y"), sep="_")
      nodes <- self$nodes
      cast.vars <- c(nodes$Lnodes,nodes$Cnodes, nodes$Anodes, nodes$Nnodes, nodes$Ynode)
      if (!missing(bslcovars)) cast.vars <- setdiff(cast.vars, bslcovars)
      odata_wide <- dcast(self$dat.sVar, formula = nodes$ID %+% " ~ " %+% nodes$tnode, value.var = cast.vars)
      return(odata_wide)
    },

    define_CVfolds = function(nfolds = 5, fold_column = "fold_id", seed = 1) {
      if (fold_column %in% names(self$dat.sVar)) {
        self$dat.sVar[, (fold_column) := NULL]
      }
      nuniqueIDs <- self$nuniqueIDs
      if (is.numeric(seed)) set.seed(seed)  #If seed is specified, set seed prior to next step
      fold_id <- sample(rep(seq(nfolds), ceiling(nuniqueIDs/nfolds)))[1:nuniqueIDs]  # Cross-validation folds (stratified folds not yet supported)
      foldsDT <- data.table("ID" = unique(self$dat.sVar[[self$nodes$IDnode]]), fold_column = fold_id)
      setnames(foldsDT, old = names(foldsDT), new = c(self$nodes$IDnode, fold_column))
      setkeyv(foldsDT, cols = self$nodes$IDnode)
      self$dat.sVar <- merge(self$dat.sVar, foldsDT, by = self$nodes$IDnode, all.x = TRUE)
      self$fold_column <- fold_column
      self$nfolds <- nfolds
      return(invisible(self))
    },

    # -----------------------------------------------------------------------------
    # (NOT USED) Create an H2OFrame and save a pointer to it as a private field
    # -----------------------------------------------------------------------------
    # load.to.H2O = function() {
    #   # if (missing(dat.sVar)) dat.sVar <- self$dat.sVar
    #   dat.sVar <- self$dat.sVar
    #   # assert_that(is.matrix(dat.sVar) | is.data.table(dat.sVar))
    #   H2O.dat.sVar <- h2o::as.h2o(dat.sVar, destination_frame = "H2O.dat.sVar")
    #   self$H2O.dat.sVar <- H2O.dat.sVar
    #   return(invisible(self))
    # },
    # -----------------------------------------------------------------------------
    # Create an H2OFrame and save a pointer to it as a private field (using faster data.table::fwrite)
    # -----------------------------------------------------------------------------
    fast.load.to.H2O = function(dat.sVar, saveH2O = TRUE, destination_frame = "H2O.dat.sVar") {
      if (missing(dat.sVar)) {
        dat.sVar <- self$dat.sVar
      }
      assert_that(is.data.table(dat.sVar))

      devDTvs <- exists("fwrite", where = "package:data.table")

      if (!devDTvs) {

        message(
"For optimal performance with h2o ML, please install the development version of data.table package.
It can be done by typing this into R terminal:
    library(devtools)
    install_github(\"Rdatatable/data.table\")")

        H2O.dat.sVar <- h2o::as.h2o(data.frame(dat.sVar), destination_frame = destination_frame)

      } else {

        types <- sapply(dat.sVar, class)
        types <- gsub("integer64", "numeric", types)
        types <- gsub("integer", "numeric", types)
        types <- gsub("double", "numeric", types)
        types <- gsub("complex", "numeric", types)
        types <- gsub("logical", "enum", types)
        types <- gsub("factor", "enum", types)
        types <- gsub("character", "string", types)
        types <- gsub("Date", "Time", types)

        # Temp file to write to:
        tmpf <- tempfile(fileext = ".csv")
        data.table::fwrite(dat.sVar, tmpf, verbose = TRUE, na = "NA_h2o")
        H2O.dat.sVar <- h2o::h2o.importFile(path = tmpf,
                                            header = TRUE,
                                            col.types = types,
                                            na.strings = rep(c("NA_h2o"), ncol(dat.sVar)),
                                            destination_frame = destination_frame)
        # replace all irregular characters to conform with destination_frame regular exprs format:
        # tmpf.dest1 <- gsub('/', 'X', tmpf, fixed = TRUE)
        # tmpf.dest2 <- gsub('.', 'X', tmpf.dest1, fixed = TRUE)
        # tmpf.dest3 <- gsub('_', 'X', tmpf.dest2, fixed = TRUE)
        # destination_frame = tmpf.dest3)
        # H2O.dat.sVar <- h2o::h2o.uploadFile(path = tmpf, parse_type = "CSV", destination_frame = "H2O.dat.sVar")
        file.remove(tmpf)

      }
      if (saveH2O) self$H2O.dat.sVar <- H2O.dat.sVar
      return(invisible(H2O.dat.sVar))
      # return(invisible(self))
    }
  ),

  active = list(
    min.t = function() { min(self$dat.sVar[[self$nodes[['tnode']]]], na.rm = TRUE) },
    max.t = function() { max(self$dat.sVar[[self$nodes[['tnode']]]], na.rm = TRUE) },
    nobs = function() { nrow(self$dat.sVar) },
    nuniqueIDs = function() { length(unique(self$dat.sVar[[self$nodes$IDnode]])) },
    nuniquets = function() { length(unique(self$dat.sVar[[self$nodes$tnode]])) },
    names.sVar = function() { colnames(self$dat.sVar) },
    ncols.sVar = function() { length(self$names.sVar) },
    dat.sVar = function(dat.sVar) {
      if (missing(dat.sVar)) {
        return(private$.mat.sVar)
      } else {
        assert_that(is.matrix(dat.sVar) | is.data.table(dat.sVar))
        private$.mat.sVar <- dat.sVar
      }
    },
   H2O.dat.sVar = function(dat.sVar) {
      if (missing(dat.sVar)) {
        return(private$.H2O.mat.sVar)
      } else {
        assert_that(is.H2OFrame(dat.sVar))
        private$.H2O.mat.sVar <- dat.sVar
      }
    },
    dat.bin.sVar = function(dat.bin.sVar) {
      if (missing(dat.bin.sVar)) {
        return(private$.mat.bin.sVar)
      } else {
        assert_that(is.matrix(dat.bin.sVar))
        private$.mat.bin.sVar <- dat.bin.sVar
      }
    },
    backup.savedGstarsDT = function() { private$.saveGstarsDT },
    emptydat.sVar = function() { private$.mat.sVar <- NULL },         # wipe out mat.sVar
    # wipe out binirized .mat.sVar:
    emptydat.bin.sVar = function() {
      private$.mat.bin.sVar <- NULL
      private$.active.bin.sVar <- NULL
    },
    noNA.Ynodevals = function(noNA.Yvals) {
      if (missing(noNA.Yvals)) return(private$.protected.YnodeVals)
      else private$.protected.YnodeVals <- noNA.Yvals
    },
    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        private$.nodes <- nodes
      }
    },
    uncensored_idx = function(uncensored_idx) {
      if (missing(uncensored_idx)) {
        return(private$.uncensored_idx)
      } else {
        assert_that(is.logical(uncensored_idx))
        private$.uncensored_idx <- uncensored_idx
      }
    },
    rule_followers_idx = function(rule_followers_idx) {
      if (missing(rule_followers_idx)) {
        return(private$.rule_followers_idx)
      } else {
        assert_that(is.logical(rule_followers_idx))
        private$.rule_followers_idx <- rule_followers_idx
      }
    },
    IPwts_by_regimen = function(IPwts) {
      if (missing(IPwts)) {
        return(private$.IPwts)
      } else {
        assert_that(is.data.table(IPwts))
        private$.IPwts <- IPwts
      }
    },
    active.bin.sVar = function() { private$.active.bin.sVar },
    ord.sVar = function() { private$.ord.sVar },
    type.sVar = function() { private$.type.sVar }
  ),

  private = list(
    .saveGstarsDT = NULL,
    .nodes = list(),              # names of the nodes in the data (Anode, Ynode, etc..)
    .protected.YnodeVals = NULL,  # Actual observed values of the binary outcome (Ynode), along with deterministic vals
    .mat.sVar = NULL,             # pointer to data.table object storing the entire dataset (including all summaries sVars)
    .H2O.mat.sVar = NULL,         # pointer to H2OFrame object that stores equivalent data to private$.mat.sVar
    .active.bin.sVar = NULL,      # Name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    .mat.bin.sVar = NULL,         # Temporary storage mat for bin indicators on currently binarized continous sVar (from private$.active.bin.sVar)
    .ord.sVar = NULL,             # Ordinal (cat) transform for continous sVar
    # sVar.object = NULL,         # DefineSummariesClass object that contains / evaluates sVar expressions
    .type.sVar = NULL,            # Named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    .uncensored_idx = NULL,       # logical vector for all observation indices that are not censored at current t
    .rule_followers_idx = NULL,   # logical vector for all observation indices that are following the current rule of interest at current t
    .IPwts = NULL,
    # Replace all missing (NA) values with a default integer (0) for matrix
    fixmiss_sVar_mat = function() {
      self$dat.sVar[gvars$misfun(self$dat.sVar)] <- gvars$misXreplace
      invisible(self)
    },
    # Replace all missing (NA) values with a default integer (0) for data.table
    fixmiss_sVar_DT = function() {
      # see http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
      dat.sVar <- self$dat.sVar
      for (j in names(dat.sVar))
        set(dat.sVar, which(gvars$misfun(dat.sVar[[j]])), j , gvars$misXreplace)
      invisible(self)
    },
    # create a vector of ordinal (categorical) vars out of cont. sVar vector:
    discretize.sVar = function(name.sVar, intervals) {
      private$.ord.sVar <- make.ordinal(x = self$get.sVar(name.sVar), intervals = intervals)
      invisible(private$.ord.sVar)
    },
    # Create a matrix of bin indicators for continuous sVar:
    binirize.cont.sVar = function(name.sVar, intervals, nbins, bin.nms) {
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = private$discretize.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms)
      invisible(self$dat.bin.sVar)
    },
    # Create a matrix of bin indicators for ordinal sVar:
    binirize.cat.sVar = function(name.sVar, levels) {
      nbins <- length(levels)
      bin.nms <- self$bin.nms.sVar(name.sVar, nbins)
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$get.sVar(name.sVar), nbins = nbins, bin.nms = bin.nms, levels = levels)
      invisible(self$dat.bin.sVar)
    },
    # Create a matrix of bin indicators for categorical CENSORING sVar (the reference category is gvars$noCENScat):
    binirize.cat.CENS.sVar = function(name.sVar, levels) {
      if (gvars$verbose) {
        message("making matrix for dummies for cat. censoring, reference category is " %+% self$noCENScat)
      }
      nbins <- length(levels)
      bin.nms <- self$bin.nms.sVar(name.sVar, nbins)
      self$dat.bin.sVar <- make.bins_mtx_cens(x.ordinal = self$get.sVar(name.sVar), nbins = nbins, bin.nms = bin.nms, levels = levels, ref.level = self$noCENScat)
      invisible(self$dat.bin.sVar)
    }
  )
)
