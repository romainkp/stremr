`[.DataStorageClass` <- function(var, indx, ...) {
  var$get.dat.sVar(indx)
}

print.DataStorageClass <- function(object){
  object$print()
}

#' Extract input data from DataStorageClass
#'
#' @param data Object of class \code{DataStorageClass} (returned by calling \code{importData} function).
#' @export
get_data = function(data) return(data$dat.sVar)

##-----------------------------------------------------------------------------
## DataStorageClass: Stores the input data;
##-----------------------------------------------------------------------------

#' Define fold ID column for cross-validation
#'
#' @param data Object of class \code{DataStorageClass} (returned by calling \code{importData} function).
#' @param nfolds The number of folds to use in V fold cross-validation.
#' @param fold_column The name for the column that will contain the fold IDs.
#' @param seed Fix the seed for random generator.
#' @export
define_CVfolds = function(data, nfolds = 5, fold_column = "fold_ID", seed = NULL) {
  data$define_CVfolds(nfolds = nfolds, fold_column = fold_column, seed = seed)
  return(data)
}

#' Removing all g model fits from DataStorage object.
#' Clean up the DataStorage object by removing all g (IPW) modeling fits.
#' This helps save up some RAM for running G-COMP / TMLE when working with large data.
#' Call this function AFTER fitting IPW (g), but before calling TMLE.
#'
#' @param OData Data storage object containing the observed data.
#' @param keep_holdout_preds Keep the holdout predictions? Needed for CV TMLE only and set to FALSE by default.
#' @export
remove_g <- function(OData, keep_holdout_preds = FALSE) {
  OData$modelfits.g0 <- NULL
  OData$modelfit.gC <- NULL
  OData$modelfit.gA <- NULL
  OData$modelfit.gN <- NULL
  if (!keep_holdout_preds) OData$g_holdout_preds <- NULL
  return(OData)
}


#' Removing not needed rows from observed data (and weights data).
#' Clean up the DataStorage object by removing all extra rows after certain time-point specified
#' by \code{tmax}. This will also have the effect of removing all the extra weights.
#' Can be useful when running TMLE with very large data, trimming the extra rows will free up some RAM.
#' @param OData Data storage object containing the observed data.
#' @param tmax Time-point value (max \code{tvals} in TMLE) after which all rows are removed.
#' @export
trim_rows_after_tmax <- function(OData, tmax = OData$max.t) {
  nodes <- OData$nodes

  keep_idx <- which(OData$dat.sVar[[nodes$tnode]] <= tmax)
  OData$dat.sVar <- OData$dat.sVar[keep_idx, ]

  if (!is.null(OData$g_preds)) {
    OData$g_preds[] <- OData$g_preds[keep_idx, ]
  }

  if (!is.null(OData$g_holdout_preds)) {
    OData$g_holdout_preds[] <- OData$g_holdout_preds[keep_idx, ]
  }

  setkeyv(OData$dat.sVar, cols = c(nodes$IDnode, nodes$tnode))
  return(OData)
}

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
#'  The pointers to this class get passed on to \code{ModelGeneric} functions: \code{$fit()},
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
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{nobs}}{...}
#'    \item{\code{ncols.sVar}}{...}
#'    \item{\code{names.sVar}}{...}
#'    \item{\code{type.sVar}}{Named list of length \code{ncol(private$.mat.sVar)} with \code{sVar} variable types: "binary"/"categor"/"contin".}
#'    \item{\code{dat.sVar}}{...}
#'    \item{\code{emptydat.sVar}}{...}
#'    \item{\code{noNA.Ynodevals}}{...}
#'    \item{\code{nodes}}{...}
#' }
#' @importFrom assertthat assert_that is.count is.flag is.string
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

    g_preds = NULL,
    g_holdout_preds = NULL,

    IPWeights_info = NULL, # A list with all the information (node names, truncation) that is needed to evaluate the IP weights.

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
    evalsubst = function(subset_vars, subset_exprs = NULL) {
      res <- rep.int(TRUE, self$nobs)
      if (!missing(subset_vars)) {
        assert_that(is.character(subset_vars))
        for (subsetvar in subset_vars) {
          # (*) find the var of interest (in self$dat.sVar or self$dat.bin.sVar), give error if not found
          sVar.vec <- self$get.outvar(var = subsetvar)
          assert_that(!is.null(sVar.vec))
          # (*) reconstruct correct expression that tests for missing values
          res <- res & (!gvars$misfun(sVar.vec))
        }
      }

      if (length(subset_exprs)==0L && !is.null(subset_exprs)) return(as.integer(subset_exprs))

      if (!is.null(subset_exprs) && !is.na(subset_exprs)) {
        if (is.logical(subset_exprs)) {
          return(which(res & subset_exprs))
        } else if (is.character(subset_exprs)) {
          ## ******************************************************
          ## data.table evaluation of the logical subset expression
          ## Note: This can be made faster by using keys in data.table on variables in eval(parse(text = subset_exprs))
          ## ******************************************************
          res.tmp <- self$dat.sVar[, eval(parse(text = subset_exprs)), by = get(self$nodes$ID)][["V1"]]
          assert_that(is.logical(res.tmp))
          return(which(res & res.tmp))
        } else if (is.integer(subset_exprs)) {
          ## The expression is already a row index, hence should be returned unchanged
          return(subset_exprs)
        }
      }
      return(which(res))
    },

    # ---------------------------------------------------------------------
    # Functions for subsetting/returning covariate design mat for ModelUnivariate Class or outcome variable
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
    # Methods for sVar types. Define the type (class) of each variable (sVar) in input data: gvars$sVartypes$bin,  gvars$sVartypes$cat or gvars$sVartypes$cont
    # --------------------------------------------------
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar), otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    # def.types.sVar = function(type.sVar = NULL) {
    #   if (is.null(type.sVar)) {
    #     private$.type.sVar <- detect.col.types(self$dat.sVar)
    #   } else {
    #     n.sVar <- length(self$names.sVar)
    #     len <- length(type.sVar)
    #     assert_that((len == n.sVar) || (len == 1L))
    #     if (len == n.sVar) { # set types for each variable
    #       assert_that(is.list(type.sVar))
    #       assert_that(all(names(type.sVar) %in% self$names.sVar))
    #     } else { # set one type for all vars
    #       assert_that(is.string(type.sVar))
    #       type.sVar <- as.list(rep(type.sVar, n.sVar))
    #       names(type.sVar) <- self$names.sVar
    #     }
    #     private$.type.sVar <- type.sVar
    #   }
    #   invisible(self)
    # },
    # set.sVar.type = function(name.sVar, new.type) { private$.type.sVar[[name.sVar]] <- new.type },
    # get.sVar.type = function(name.sVar) { if (missing(name.sVar)) { private$.type.sVar } else { private$.type.sVar[[name.sVar]] } },
    # is.sVar.bin = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$bin },
    # is.sVar.cat = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$cat },
    # is.sVar.cont = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$cont },
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
      nodes <- nodes[!is.na(nodes)]
      for (node in nodes) CheckVarNameExists(self$dat.sVar, node)
      private$.saveGstarsDT <- self$dat.sVar[, nodes, with = FALSE]
      return(invisible(self))
    },

    restoreNodes = function(nodes) {
      nodes <- nodes[!is.null(nodes)]
      nodes <- nodes[!is.na(nodes)]
      for (node in nodes) CheckVarNameExists(self$dat.sVar, node)
      if (is.null(private$.saveGstarsDT)) stop("Nodes in dat.sVar cannot be restored, private$.saveGstarsDT is null!")
      self$dat.sVar[, (nodes) := private$.saveGstarsDT[, nodes, with = FALSE]]
      return(invisible(self))
    },

    ## Modify the values in node nodes_to_repl in self$dat.sVar with values from source_for_repl using only the IDs in subset_idx
    ## This is done by reference (modifying the input data.table)
    replaceNodesVals = function(subset_idx, nodes_to_repl = intervened_NODE, source_for_repl = NodeNames) {
      for (node_idx in seq_along(nodes_to_repl)) {
        if (length(subset_idx) > 0) {
          source_node <- self$dat.sVar[subset_idx, (source_for_repl[node_idx]), with = FALSE][[source_for_repl[node_idx]]]
          self$dat.sVar[subset_idx, (nodes_to_repl[node_idx]) := source_node]
        }
      }
      return(invisible(self))
    },

    ## Rescale the node values (multply by delta) in a data.table by reference
    rescaleNodes = function(subset_idx, nodes_to_rescale, delta) {
      assert_that(is.numeric(delta) || is.integer(delta))
      for (node_idx in seq_along(nodes_to_rescale)) {
        if (length(subset_idx) > 0) {
          self$dat.sVar[subset_idx, (nodes_to_rescale[node_idx]) := get(nodes_to_rescale[node_idx])*delta]
        }
      }
      return(invisible(self))
    },

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

    ## Rule followers over time for given exposure node and the counterfactual intervention node
    ## ******* Must follow the rule over all time-points up to and including t *******
    ## This is used for stratifying GCOMP by the entire history of following teh rule.
    ## For stochastic intervention nodes, the person follows the rule when the joint probability of g^*(t) over all t is  > 0.
    ## Returns a logical vector of indicators of length nrow(data)
    eval_follow_rule = function(NodeName, gstar.NodeName) {
      dat.rulefollow <- self$dat.sVar[, c(self$nodes$IDnode, NodeName, gstar.NodeName), with = FALSE]
      setkeyv(dat.rulefollow, cols = self$nodes$IDnode)
      ## Likelihood of N(t) under counterfactual g^*. This is done so that stochastic intervention nodes are handled appropriately.
      dat.rulefollow[, "follow.byt" := get(gstar.NodeName)^get(NodeName) * (1L - get(gstar.NodeName))^(1 - get(NodeName))]
      ## rule follower whenever the probability of observing the current history \bar{(O(t))} is > 0
      followers <- dat.rulefollow[, "followers" := as.logical(cumprod(follow.byt) > 0), by = eval(self$nodes$IDnode)][["followers"]]
      return(followers)
    },

    ## Time-specific rule followers for given exposure node and the counterfactual intervention node
    ## **** Only considers the equality of the exposure and intervention nodes at the current t, regardless of the past history ***
    ## This is used for stratifying GCOMP by current time-point only.
    ## For stochastic intervention nodes, the person follows the rule when the probability of g^* at t is > 0.
    ## Returns a logical vector of indicators of length nrow(data)
    eval_follow_rule_each_t = function(NodeName, gstar.NodeName) {
      dat.rulefollow <- self$dat.sVar[, c(self$nodes$IDnode, NodeName, gstar.NodeName), with = FALSE]
      setkeyv(dat.rulefollow, cols = self$nodes$IDnode)
      ## likelihood of N(t) under counterfactual g^*.
      ## Rule follower at t whenever the probability of observing current gstar(t) is > 0.
      follower_by_t <- dat.rulefollow[,
                          "follow.byt" := as.logical(get(gstar.NodeName)^get(NodeName) * (1L - get(gstar.NodeName))^(1 - get(NodeName)) > 0)][[
                          "follow.byt"]]
      return(follower_by_t)
    },

    ## logical vector of length nrow(data), indicating if the person has not been censored (yet) at current t
    eval_uncensored = function() {
      if (!is.null(self$nodes$Cnodes)) {
        return(self$dat.sVar[, list(uncensored = as.logical(rowSums(.SD, na.rm = TRUE) == eval(self$noCENScat))), .SDcols = self$nodes$Cnodes][["uncensored"]])
      } else {
        return(rep.int(TRUE, nrow(self$dat.sVar)))
      }
    },

    ## integer vector of rows in data for observations (over all time-points) who have not been censored yet at each t
    eval_uncensored_idx = function() {
      return(self$dat.sVar[,
              list(uncensored_idx = which(as.logical(rowSums(.SD, na.rm = TRUE) == eval(self$noCENScat)))),
              .SDcols = self$nodes$Cnodes][[
                "uncensored_idx"]])
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
      return(invisible(self))
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

    define_CVfolds = function(nfolds = 5, fold_column = "fold_ID", seed = NULL) {

      ## Means that fold column is already defined in the input data, just copy the information and perform few checks
      if (missing(nfolds)) {
        if (missing(fold_column)) stop("fold_column must be specified when nfolds is missing")
        if (!fold_column %in% names(self$dat.sVar)) stop("fold_column could not be located in the input data")
        if (!is.integerish(self$dat.sVar[[fold_column]]) && !(is.factor(self$dat.sVar[[fold_column]])))
          stop("'fold_column must be either an integer or a factor")
        nfolds <- length(unique(self$dat.sVar[[fold_column]]))

        self$fold_column <- fold_column
        ## evaluate the number of unique folds:
        self$nfolds <- nfolds
        return(invisible(self))

      } else {

        if (fold_column %in% names(self$dat.sVar)) self$dat.sVar[, (fold_column) := NULL]
        nuniqueIDs <- self$nuniqueIDs

        if (is.numeric(seed)) set.seed(seed)  #If seed is specified, set seed prior to next step

        fold_IDs <- sprintf("%02d", seq(nfolds))
        fold_id <- as.factor(sample(rep(fold_IDs, ceiling(nuniqueIDs/nfolds)))[1:nuniqueIDs])  # Cross-validation folds
        # fold_id <- sample(rep(seq(nfolds), ceiling(nuniqueIDs/nfolds)))[1:nuniqueIDs]  # Cross-validation folds (stratified folds not yet supported)

        foldsDT <- data.table("ID" = unique(self$dat.sVar[[self$nodes$IDnode]]), fold_column = fold_id)
        setnames(foldsDT, old = names(foldsDT), new = c(self$nodes$IDnode, fold_column))
        setkeyv(foldsDT, cols = self$nodes$IDnode)
        self$dat.sVar <- merge(self$dat.sVar, foldsDT, by = self$nodes$IDnode, all.x = TRUE)
        self$fold_column <- fold_column
        self$nfolds <- nfolds

        if (is.numeric(seed)) set.seed(NULL)  #If seed was specified, reset it to NULL
        return(invisible(self))
      }
    },

    make_origami_fold_from_column = function(subset_idx) {
      fold_column <- self$fold_column
      if (is.null(fold_column)) stop("must define the column with validation folds")
      # n <- nrow(data)
      n <- length(subset_idx)
      folds <- self$dat.sVar[subset_idx, ][[fold_column]]
      k <- length(unique(folds))
      idx <- seq_len(n)
      fold_idx <- split(idx, folds)
      fold <- function(v, test) {
        origami::make_fold(v, setdiff(idx, test), test)
      }
      purrr::map2((1:k), fold_idx, fold)
    },

    type.sVar = function(var_names, pcontinuous = 0.05) {
      types_list <- lapply(var_names, function(var) {
        x <- self$get.sVar(var)
        x <- x[!is.na(x)] 
        var_obj <- sl3::variable_type(x = x, pcontinuous = pcontinuous)
        var_obj$type
      })
      return(types_list)
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
    backup.savedGstarsDT = function() { private$.saveGstarsDT },
    emptydat.sVar = function() { private$.mat.sVar <- NULL },         # wipe out mat.sVar
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
    uncensored = function(uncensored) {
      if (missing(uncensored)) {
        return(private$.uncensored)
      } else {
        assert_that(is.logical(uncensored))
        private$.uncensored <- uncensored
      }
    },
    follow_rule = function(follow_rule) {
      if (missing(follow_rule)) {
        return(private$.follow_rule)
      } else {
        assert_that(is.logical(follow_rule))
        private$.follow_rule <- follow_rule
      }
    },
    IPwts_by_regimen = function(IPwts) {
      if (missing(IPwts)) {
        return(private$.IPwts)
      } else {
        assert_that(is.data.table(IPwts))
        private$.IPwts <- IPwts
      }
    }
  ),

  private = list(
    .saveGstarsDT = NULL,
    .nodes = list(),              # names of the nodes in the data (Anode, Ynode, etc..)
    .protected.YnodeVals = NULL,  # Actual observed values of the binary outcome (Ynode), along with deterministic vals
    .mat.sVar = NULL,             # pointer to data.table object storing the entire dataset (including all summaries sVars)
    .active.bin.sVar = NULL,      # Name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    .mat.bin.sVar = NULL,         # Temporary storage mat for bin indicators on currently binarized continous sVar (from private$.active.bin.sVar)
    .ord.sVar = NULL,             # Ordinal (cat) transform for continous sVar
    # .type.sVar = NULL,            # Named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    .uncensored = NULL,       # logical vector for all observation indices that are not censored at current t
    .follow_rule = NULL,   # logical vector for all observation indices that are following the current rule of interest at current t
    .IPwts = NULL
    # Replace all missing (NA) values with a default integer (0) for matrix
    # fixmiss_sVar_mat = function() {
    #   self$dat.sVar[gvars$misfun(self$dat.sVar)] <- gvars$misXreplace
    #   invisible(self)
    # },
    # Replace all missing (NA) values with a default integer (0) for data.table
    # fixmiss_sVar_DT = function() {
    #   # see http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
    #   dat.sVar <- self$dat.sVar
    #   for (j in names(dat.sVar))
    #     set(dat.sVar, which(gvars$misfun(dat.sVar[[j]])), j , gvars$misXreplace)
    #   invisible(self)
    # }
  )
)
