## ------------------------------------------------------------------------------------------
## TO DO:
## ------------------------------------------------------------------------------------------
## *) The definition of Qperiods below needs to be based on actual periods observed in the data.
## *) Stochastic interventions ***
##    intervened_TRT/intervened_MONITOR are vectors of counterfactual probs -> need to allow each to be multivariate (>1 cols)
## *) Accessing correct ModelQlearn ***
##     Need to come up with a good way of accessing correct terminal ModelQlearn to get final Q-predictions in private$probAeqa.
##     These predictions are for the final n observations that were used for evaluating E[Y^a].
##     However, if we do stratify=TRUE and use a stack version of g-comp with different strata all in one stacked database it will be difficult
##     to pull the right ModelQlearn.
##     Alternative is to ignore that completely and go directly for the observed data (by looking up the right column/row combos).
## *) Pooling across regimens ***
##     Since new regimen results in new Qkplus1 and hence new outcome -> requires a separate regression for each regimen:
##     => Can either pool all Q.regimens at once (fit one model for repeated observations, one for each rule, smoothing over rules).
##     => Can use the same stacked dataset of regime-followers, but with a separate stratification for each regimen.
## ------------------------------------------------------------------------------------------

## ------------------------------------------------------------------------------------------
## Sequential (Recursive) G-COMP Algorithm:
## ------------------------------------------------------------------------------------------
## 1. At each t iteration:
##     Fitting:
##        Outcome is always Qkplus1 (either by regimen or not). Fitting subset excludes all obs censored at t.
##     Prediction:
##        Swap the entire column(s) A[t] with A^*[t] by renaming them in OData$dat.sVar (might do this for more than one regimen in the future, if pooling Q's)
##        Subset all observation by t (including all censored and non-followers).
##        PredictP1() for everybody in the subset, then save the prediction in rows Qkplus1[t] and rows Qkplus1[t-1] for all obs that were used in prediction.
##        Results in correct algorithm even when stratifyQ_by_rule=TRUE: next iteration (t-1) will do fit only based on obs that followed the rule at t-1.
## 2. At next iteration t-1:
##     Fitting:
##        Use the predictions in Qkplus1[t-1] as new outcomes and repeat until reached minimum t.

## ------------------------------------------------------------------------------------------
## Stochastic interventions:
## ------------------------------------------------------------------------------------------
## When the intervention node values are <1 and >0, those are interpreted as probabilities,
## i.e., stochastic interventions.
## In more detail, the column 0 < intervened_TRT < 1 defines the counterfactual probability that P(TRT[t]^*=1)=intervened_TRT[t].
## Those are dealt with in ModelQlearn class by directly integrating, i.e.,
##  evaluating a weighted sum of predicted Q's with weights given by the probabilities in g.star columns
## (on A and N).

## ------------------------------------------------------------------------------------------
##  *** rule_followers:
## ------------------------------------------------------------------------------------------
##     Once you go off the treatment first time, this is it, the person is censored for the rest of the follow-up.
##     Rule followers are evaluated automatically by comparing (intervened_TRT and TRT and CENS) and (intervened_MONITOR and MONITOR and CENS).
##     Rule followers are: (intervened_TRT = 1) & (TRT == 1) or (intervened_TRT = 0) & (TRT == 0) or (intervened_TRT > 0 & intervened_TRT < 0) & (Not Censored).
##     Exactly the same logic is also applied to (intervened_MONITOR & MONITOR) when these are specified.
##     NOT TESTED: If either TRT or MONITOR is multivariate (>1 col), this logic needs to be applied FOR EACH COLUMN of TRT/MONITOR.
## ------------------------------------------------------------------------------------------
## *** Pooling & performing only one TMLE update across all Q.k ***
## ------------------------------------------------------------------------------------------
##     In general, we can pool all Q.kplus (across all t's) and do a single update on all of them
##     Would make this TMLE iterative, but would not require refitting on the initial Q's


## ------------------------------------------------------------------------------------------------------------------------
## When useonly_t_TRT or useonly_t_MONITOR is specified, need to set nodes to their observed values, rather than the counterfactual values
## ------------------------------------------------------------------------------------------------------------------------
## Do it separately for gstar_TRT & gstar_MONITOR
## Loop over each node in gstar_TRT / gstar_MONITOR
## Do it only once for all observations inside main TMLE call
## Back-up a copy of all gstar nodes first, the original copy is then restored when finished running
## The observations which get swapped with g0 values are defined by:
## subset_idx <- OData$evalsubst(subset_exprs = useonly_t_NODE)
## probability of P(A^*(t)=n(t)) or P(N^*(t)=n(t)) under counterfactual A^*(t) or N^*(t) and observed a(t) or n(t)
## Example call: defineNodeGstarGCOMP(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, stratifyQ_by_rule)
defineNodeGstarGCOMP <- function(OData, intervened_NODE, NodeNames, useonly_t_NODE, stratifyQ_by_rule, stratify_by_last) {
  # if intervened_NODE returns more than one rule-column, evaluate g^* for each and the multiply to get a single joint (for each time point)
  if (!is.null(intervened_NODE) && !is.na(intervened_NODE)) {
    gstar.NODEs <- intervened_NODE
    for (intervened_NODE_col in intervened_NODE) CheckVarNameExists(OData$dat.sVar, intervened_NODE_col)
    assert_that(length(intervened_NODE) == length(NodeNames))

    ## ------------------------------------------------------------------------------------------
    ## Modify the observed input intervened_NODE in OData$dat.sVar with values from NodeNames for subset_idx:
    ## ------------------------------------------------------------------------------------------
    subset_idx <- OData$evalsubst(subset_exprs = useonly_t_NODE)
    not_subset_idx <- setdiff(1:nrow(OData$dat.sVar), subset_idx)
    # OData$replaceNodesVals(!subset_idx, nodes_to_repl = intervened_NODE, source_for_repl = NodeNames)
    OData$replaceNodesVals(not_subset_idx, nodes_to_repl = intervened_NODE, source_for_repl = NodeNames)
    ## ------------------------------------------------------------------------------------------

    ## ------------------------------------------------------------------------------------------
    ## update rule followers for trt if doing stratified G-COMP:
    ## Note this will define rule followers based on REPLACED intervened_NODE in dat.sVar (i.e., modified n^*(t) under N.D.E.)
    ## FOR NDE BASED TMLE THE DEFINITION OF RULE-FOLLOWERS CHANGES ACCORDINGLY based on modified n^*(t) and a^*(t)
    # ------------------------------------------------------------------------------------------
    if (stratifyQ_by_rule) {
      if (!stratify_by_last) {
        follow_rule <- OData$eval_follow_rule(NodeName = NodeNames, gstar.NodeName = intervened_NODE)
      } else {
        follow_rule <- OData$eval_follow_rule_each_t(NodeName = NodeNames, gstar.NodeName = intervened_NODE)
      }
      OData$follow_rule <- follow_rule & OData$follow_rule & OData$uncensored
    }
  } else {
    # use the actual (observed) node names under g0:
    gstar.NODEs <- NodeNames
  }
  return(gstar.NODEs)
}

# ---------------------------------------------------------------------------------------
#' Iterative TMLE wrapper for \code{fit_GCOMP}
#'
#' Calls \code{fit_GCOMP} with argument \code{iterTMLE = TRUE}.
#' @param ... Arguments that will be passed down to the underlying function \code{fit_GCOMP}
#' @return \code{data.table} with survival by time for sequential GCOMP and iterative TMLE
#' @seealso \code{\link{fit_GCOMP}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fit_iterTMLE <- function(...) {
  fit_GCOMP(TMLE = FALSE, iterTMLE = TRUE, ...)
}

# ---------------------------------------------------------------------------------------
#' TMLE wrapper for \code{fit_GCOMP}
#'
#' Calls \code{fit_GCOMP} with argument \code{TMLE = TRUE}.
#' @param ... Arguments that will be passed down to the underlying function \code{fit_GCOMP}
#' @return \code{data.table} with TMLE survival by time
#' @seealso \code{\link{fit_GCOMP}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fit_TMLE <- function(...) {
  fit_GCOMP(TMLE = TRUE, ...)
}

# ---------------------------------------------------------------------------------------
#' CV-TMLE wrapper for \code{fit_GCOMP}
#'
#' Calls \code{fit_GCOMP} with arguments \code{TMLE = TRUE} and \code{CVTMLE = TRUE}.
#' @param ... Arguments that will be passed down to the underlying function \code{fit_GCOMP}
#' @return \code{data.table} with TMLE survival by time
#' @seealso \code{\link{fit_GCOMP}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fit_CVTMLE <- function(...) {
  fit_GCOMP(TMLE = TRUE, CVTMLE = TRUE, ...)
}


# ---------------------------------------------------------------------------------------
#' Fit sequential GCOMP and TMLE for survival
#'
#' Interventions on up to 3 nodes are allowed: \code{CENS}, \code{TRT} and \code{MONITOR}.
#' TMLE adjustment will be based on the inverse of the propensity score fits for the observed likelihood (g0.C, g0.A, g0.N),
#' multiplied by the indicator of not being censored and the probability of each intervention in \code{intervened_TRT} and \code{intervened_MONITOR}.
#' Requires column name(s) that specify the counterfactual node values or the counterfactual probabilities of each node being 1 (for stochastic interventions).
#' @param OData Input data object created by \code{importData} function.
#' @param tvals Vector of time-points in the data for which the survival function (and risk) should be estimated
#' @param Qforms Regression formulas, one formula per Q. Only main-terms are allowed.
#' @param Qstratify Placeholder for future user-defined model stratification for fitting Qs (CURRENTLY NOT FUNCTIONAL, WILL RESULT IN ERROR).
#' @param intervened_TRT Column name in the input data with the probabilities (or indicators) of counterfactual treatment nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on treatment node(s).
#' @param intervened_MONITOR Column name in the input data with probabilities (or indicators) of counterfactual monitoring nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on the monitoring node(s).
#' @param useonly_t_TRT Use for intervening only on some subset of observation and time-specific treatment nodes.
#' Should be a character string with a logical expression that defines the subset of intervention observations.
#' For example, using \code{TRT==0} will intervene only at observations with the value of \code{TRT} being equal to zero.
#' The expression can contain any variable name that was defined in the input dataset.
#' Leave as \code{NULL} when intervening on all observations/time-points.
#' @param useonly_t_MONITOR Same as \code{useonly_t_TRT}, but for monitoring nodes.
#' @param rule_name Optional name for the treatment/monitoring regimen.
#' @param stratifyQ_by_rule Set to \code{TRUE} for stratifying the fit of Q (the outcome model) by rule-followers only.
#' There are two ways to do this stratification. The first option is to use \code{stratify_by_last=TRUE}  (default),
#' which would fit the outcome model only among the observations that were receiving their supposed
#' counterfactual treatment at the current time-point (ignoring the past history of treatments leading up to time-point t).
#' The second option is to set \code{stratify_by_last=FALSE} in which case the outcome model will be fit only
#' among the observations who followed their counterfactual treatment regimen throughout the entire treatment history up to
#' current time-point t (rule followers). For the latter option, the observation would be considered a non-follower if
#' the person's treatment did not match their supposed counterfactual treatment at any time-point up to and including current
#' time-point t.
#' @param stratify_by_last Only used when \code{stratifyQ_by_rule} is \code{TRUE}.
#' Set to \code{TRUE} for stratification by last time-point, set to \code{FALSE} for stratification by all time-points (rule-followers).
#' See \code{stratifyQ_by_rule} for more details.
#' @param TMLE Set to \code{TRUE} to run the usual longitudinal TMLE algorithm (with a separate TMLE update of Q for every sequential regression).
#' @param iterTMLE Set to \code{TRUE} to run the iterative univariate TMLE instead of the usual longitudinal TMLE.
#' When set to \code{TRUE} this will also provide the standard sequential Gcomp as party of the output.
#' @param CVTMLE Set to \code{TRUE} to run the CV-TMLE algorithm instead of the usual TMLE algorithm.
#' Must set either \code{TMLE}=\code{TRUE} or \code{iterTMLE}=\code{TRUE} for this argument to have any effect.
#' @param byfold_Q (ADVANCED USE) Fit iterative means (Q parameter) using "by-fold" (aka "fold-specific" or "split-specific") cross-validation approach.
#' Only works with \code{fit_method}=\code{"origamiSL"}.
#' @param IPWeights (Optional) result of calling function \code{getIPWeights} for running TMLE (evaluated automatically when missing)
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @param models Optional parameters specifying the models for fitting the iterative (sequential) G-Computation formula.
#' Must be an object of class \code{ModelStack} specified with \code{gridisl::defModel} function.
#' @param weights Optional \code{data.table} with additional observation-time-specific weights.  Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param max_iter For iterative TMLE only: Integer, set to maximum number of iterations for iterative TMLE algorithm.
#' @param adapt_stop For iterative TMLE only: Choose between two stopping criteria for iterative TMLE, default is \code{TRUE},
#' which will stop the iterative TMLE algorithm in an adaptive way. Specifically, the iterations will stop when the mean estimate
#' of the efficient influence curve is less than or equal to 1 / (\code{adapt_stop_factor}*sqrt(\code{N})), where
#' N is the total number of unique subjects in data and \code{adapt_stop_factor} is set to 10 by default.
#' When \code{TRUE}, the argument \code{tol_eps} is ignored and TMLE stops when either \code{max_iter} has been reached or this criteria has been satisfied.
#' When \code{FALSE}, the stopping criteria is determined by values of \code{max_iter} and \code{tol_eps}.
#' @param adapt_stop_factor For iterative TMLE only: The adaptive factor to choose the stopping criteria for iterative TMLE when
#' \code{adapt_stop} is set to \code{TRUE}. Default is 10.
#' TMLE will keep iterative until
#' the mean estimate of the efficient influence curve is less than 1 / (\code{adapt_stop_factor}*sqrt(\code{N})) or when the number of iterations is \code{max_iter}.
#' @param tol_eps For iterative TMLE only: Numeric error tolerance for the iterative TMLE update.
#' The iterative TMLE algorithm will stop when the absolute value of the TMLE intercept update is below \code{tol_eps}
#' @param parallel Set to \code{TRUE} to run the sequential G-COMP or TMLE in parallel (uses \code{foreach} with \code{dopar} and
#' requires a previously defined parallel back-end cluster)
#' @param estimator Specify the default estimator to use for fitting the iterative g-computation formula.
#' Should be a character string in the format 'Package__Algorithm'.
#' See \code{stremrOptions("estimator", showvals = TRUE)} for a range of possible values.
#' This argument is ignored when the fitting procedures are already defined via the argument \code{models}.
#' @param fit_method Model selection approach. Can be either \code{"none"} - no model selection or
#' \code{"cv"} - discrete Super Learner using V fold cross-validation that selects the best model according to lowest cross-validated MSE (must specify the column name that contains the fold IDs) or
#' \code{"origamiSL"} - continuous Super Learner that uses the \code{origami} R package to select the
#' convex combination of the model predictions (aka model stacking).
#  \code{"holdout"} - model selection by splitting the data into training and validation samples according to lowest validation sample MSE (must specify the column of \code{TRUE} / \code{FALSE} indicators,
#  where \code{TRUE} indicates that this row will be selected as part of the model validation sample).
#' @param fold_column The column name in the input data (ordered factor) that contains the fold IDs to be used as part of the validation sample.
#' Use the provided function \code{\link{define_CVfolds}} to
#' define such folds or define the folds using your own method.
#' @param return_wts Applies only when \code{TMLE = TRUE}.
#' Return the data.table with subject-specific IP weights as part of the output.
#' Note: for large datasets setting this to \code{TRUE} may lead to extremely large object sizes!
#' @param return_fW Return the \code{gridisl} model object from the very last Q regression.
#' Can be used for obtaining subject-specific predictions of the counterfactual functional E(Y_{d}|W_i).
#' @param reg_Q (ADVANCED USE ONLY) Directly specify the Q regressions, separately for each time-point.
#' @param type_intervened_TRT (ADVANCED FEATURE) TBD
#' @param type_intervened_MONITOR (ADVANCED FEATURE) TBD
#' @param maxpY Maximum probability that the cumulative incidence of the outcome Y(t) is equal to 1.
#' Useful for upper-bounding the rare-outcomes.
#' @param TMLE_updater Function for performing the TMLE update. Default is the TMLE updater based on speedglm (called \code{"TMLE.updater.speedglm"}).
#' Other possible options include \code{"TMLE.updater.glm"}, \code{"linear.TMLE.updater.speedglm"} and \code{"iTMLE.updater.xgb"}.
#' @param verbose Set to \code{TRUE} to print auxiliary messages during model fitting.
#' @param ... When \code{models} arguments is NOT specified, these additional arguments will be passed on directly to all \code{GridSL}
#' modeling functions that are called from this routine,
#' e.g., \code{family = "binomial"} can be used to specify the model family.
#' Note that all such arguments must be named.
#' @return An output list containing the \code{data.table} with survival estimates over time saved as \code{"estimates"}.
#' @seealso \code{\link{stremr-package}} for the general overview of the package.
#' @example tests/examples/2_building_blocks_example.R
#' @export
# proposed name change:
fit_GCOMP <- function(OData,
                        tvals,
                        Qforms,
                        intervened_TRT = NULL,
                        intervened_MONITOR = NULL,
                        rule_name = paste0(c(intervened_TRT, intervened_MONITOR), collapse = ""),
                        models = NULL,
                        estimator = stremrOptions("estimator"),
                        fit_method = stremrOptions("fit_method"),
                        fold_column = stremrOptions("fold_column"),
                        TMLE = FALSE,
                        stratifyQ_by_rule = FALSE,
                        stratify_by_last = TRUE,
                        Qstratify = NULL,
                        useonly_t_TRT = NULL,
                        useonly_t_MONITOR = NULL,
                        iterTMLE = FALSE,
                        CVTMLE = FALSE,
                        byfold_Q = FALSE,
                        IPWeights = NULL,
                        trunc_weights = 10^6,
                        weights = NULL,
                        max_iter = 15,
                        adapt_stop = TRUE,
                        adapt_stop_factor = 10,
                        tol_eps = 0.001,
                        parallel = FALSE,
                        return_wts = FALSE,
                        return_fW = FALSE,
                        reg_Q = NULL,
                        type_intervened_TRT = NULL,
                        type_intervened_MONITOR = NULL,
                        maxpY = 1.0,
                        TMLE_updater = "TMLE.updater.speedglm",
                        verbose = getOption("stremr.verbose"), ...) {

  # cat("Calling fit_GCOMP:\n")
  # cat("intervened_TRT: ", intervened_TRT, "\n")
  # cat("stratifyQ_by_rule: ", stratifyQ_by_rule, "\n")
  # cat("stratify_by_last: ", stratify_by_last, "\n")
  # cat("trunc_weights: ", trunc_weights, "\n")

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  if (!is.null(models)) {
    assert_that(is.ModelStack(models) || is(models, "Lrnr_base"))
  }

  assert_that(is.logical(adapt_stop))

  if (TMLE & iterTMLE) stop("Either 'TMLE' or 'iterTMLE' must be set to FALSE. Cannot estimate both within a single algorithm run.")

  if (byfold_Q && !(fit_method %in% "origamiSL"))
    stop("Running byfold_Q=TRUE requires seting fit_method='origamiSL'")

  if (missing(rule_name)) rule_name <- paste0(c(intervened_TRT,intervened_MONITOR), collapse = "")

  added_rule_col <- FALSE
  if (!"rule.name" %in% names(OData$dat.sVar)) {
    added_rule_col <- TRUE
    rule_name <- paste0(c(intervened_TRT,intervened_MONITOR), collapse = "")
    OData$dat.sVar[, "rule.name" := rule_name]
  } else {
    rule_name <- unique(OData$dat.sVar[["rule.name"]])
  }

  # ------------------------------------------------------------------------------------------------
  # **** Evaluate the uncensored and initialize rule followers (everybody is a follower by default)
  # ------------------------------------------------------------------------------------------------
  OData$uncensored <- OData$eval_uncensored()
  OData$follow_rule <- rep.int(TRUE, nrow(OData$dat.sVar)) # (everybody is a follower by default)

  sVar.exprs <- capture.exprs(...)
  models_control <- c(list(models = models), list(reg_Q = reg_Q), opt_params = list(sVar.exprs))
  models_control[["estimator"]] <- estimator[1L]
  models_control[["fit_method"]] <- fit_method[1L]
  models_control[["fold_column"]] <- fold_column

  # ------------------------------------------------------------------------------------------------
  # **** Add weights if TMLE=TRUE and if weights were defined
  # NOTE: This needs to be done only once if evaluating survival over several time-points
  # ------------------------------------------------------------------------------------------------
  if (TMLE || iterTMLE) {
    if (is.null(IPWeights)) {
      if (gvars$verbose) message("...evaluating IPWeights for TMLE...")
      IPWeights <- getIPWeights(OData,
                                intervened_TRT,
                                intervened_MONITOR,
                                useonly_t_TRT,
                                useonly_t_MONITOR,
                                rule_name,
                                holdout = CVTMLE,
                                eval_stabP = FALSE,
                                type_intervened_TRT = type_intervened_TRT,
                                type_intervened_MONITOR = type_intervened_MONITOR)
      if (trunc_weights < Inf) IPWeights[eval(as.name("cum.IPAW")) > trunc_weights, ("cum.IPAW") := trunc_weights]
    } else {
      getIPWeights_fun_call <- attributes(IPWeights)[['getIPWeights_fun_call']]
      if (gvars$verbose) message("applying user-specified IPWeights obtained with a call: \n" %+% deparse(getIPWeights_fun_call)[[1]])
      assert_that(all.equal(attributes(IPWeights)[['intervened_TRT']], intervened_TRT))
      assert_that(all.equal(attributes(IPWeights)[['intervened_MONITOR']], intervened_MONITOR))
    }
    assert_that(is.data.table(IPWeights))
    assert_that("cum.IPAW" %in% names(IPWeights))
    OData$IPwts_by_regimen <- IPWeights

    if (!is.null(weights)) stop("optional argument 'weights' is not implemented for TMLE or GCOMP")
    ## Add additional observation-specific weights to the cumulative weights:
    # IPWeights <- process_opt_wts(IPWeights, weights, nodes, adjust_outcome = FALSE)
  }

  if (missing(tvals)) stop("must specify survival 'tvals' of interest (time period values from column " %+% nodes$tnode %+% ")")

  # ------------------------------------------------------------------------------------------
  # Create a back-up of the observed input gstar nodes (created by user in input data):
  # Will add new columns (that were not backed up yet) TO SAME backup data.table
  # ------------------------------------------------------------------------------------------
  OData$backupNodes(c(intervened_TRT,intervened_MONITOR))

  # ------------------------------------------------------------------------------------------------
  # Define the intervention nodes
  # Modify the observed input intervened_NODE in OData$dat.sVar with values from NodeNames for subset_idx
  # ------------------------------------------------------------------------------------------------
  gstar.A <- defineNodeGstarGCOMP(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, stratifyQ_by_rule, stratify_by_last)
  gstar.N <- defineNodeGstarGCOMP(OData, intervened_MONITOR, nodes$Nnodes, useonly_t_MONITOR, stratifyQ_by_rule, stratify_by_last)
  interventionNodes.g0 <- c(nodes$Anodes, nodes$Nnodes)
  interventionNodes.gstar <- c(gstar.A, gstar.N)

  # When the same node names belongs to both g0 and gstar it doesn't need to be intervened upon, so exclude
  common_names <- intersect(interventionNodes.g0, interventionNodes.gstar)
  interventionNodes.g0 <- interventionNodes.g0[!interventionNodes.g0 %in% common_names]
  interventionNodes.gstar <- interventionNodes.gstar[!interventionNodes.gstar %in% common_names]

  OData$interventionNodes.g0 <- interventionNodes.g0
  OData$interventionNodes.gstar <- interventionNodes.gstar

  ## ------------------------------------------------------------------------------------------------
  ## RUN GCOMP OR TMLE FOR SEVERAL TIME-POINTS EITHER IN PARALLEL OR SEQUENTIALLY
  ## For est of S(t) over vector of ts, estimate for highest t first going down to smallest t
  ## This is a more efficient when parallelizing, since larger t implies more model runs & longer run time
  ## ------------------------------------------------------------------------------------------------
  est_name <- ifelse(TMLE, "TMLE", ifelse(iterTMLE, "GCOMP & iterTMLE", "GCOMP"))

  ## using delayed (works with multisession, errors with multcore)
  subtasks <- lapply(rev(seq_along(tvals)), function(t_idx) {
    t_period <- tvals[t_idx]
    delayed_Q <- delayed_fun(fit_GCOMP_onet)(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule,
                   TMLE = TMLE, iterTMLE = iterTMLE, CVTMLE = CVTMLE, byfold_Q = byfold_Q,
                   models = models_control, max_iter = max_iter, adapt_stop = adapt_stop,
                   adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps,
                   return_fW = return_fW, maxpY = maxpY, TMLE_updater = TMLE_updater, verbose = verbose)
    delayed_Q$expect_error <- TRUE
    return(delayed_Q)
  })
  bundle_Q <- bundle_delayed(subtasks)

  tmle.run.res <- try(
    res_byt <- bundle_Q$compute(verbose = TRUE)
    # res_byt <- bundle_Q$compute(nworkers = 2, verbose = TRUE)
  )

  ## using native future interface (error: result is too long a vector):
  # f <- list()
  # for (ii in seq_along(tvals)) {
  #   t_period <- tvals[ii]
  #   f[[ii]] <- future({
  #     fit_GCOMP_onet(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule,
  #                    TMLE = TMLE, iterTMLE = iterTMLE, CVTMLE = CVTMLE, byfold_Q = byfold_Q,
  #                    models = models_control, max_iter = max_iter, adapt_stop = adapt_stop,
  #                    adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps,
  #                    return_fW = return_fW, maxpY = maxpY, TMLE_updater = TMLE_updater, verbose = verbose)
  #   })
  # }
  # #
  # res <- f
  # res <- lapply(f, FUN = value)

  ## previous version with foreach() paradigm:
  # tmle.run.res <- try(
    # res_byt <- foreach::foreach(t_idx = rev(seq_along(tvals)), .options.multicore = mcoptions) %dopar% {
    #   t_period <- tvals[t_idx]
    #   res <- fit_GCOMP_onet(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule,
    #                           TMLE = TMLE, iterTMLE = iterTMLE, CVTMLE = CVTMLE, byfold_Q = byfold_Q,
    #                           models = models_control, max_iter = max_iter, adapt_stop = adapt_stop,
    #                           adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps,
    #                           return_fW = return_fW, maxpY = maxpY, TMLE_updater = TMLE_updater, verbose = verbose)
    #   return(res)
    # if (parallel) {
    #   mcoptions <- list(preschedule = FALSE)
    #   '%dopar%' <- foreach::'%dopar%'
      # res_byt <- foreach::foreach(t_idx = rev(seq_along(tvals)), .options.multicore = mcoptions) %dopar% {
      #   t_period <- tvals[t_idx]
      #   res <- fit_GCOMP_onet(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule,
      #                           TMLE = TMLE, iterTMLE = iterTMLE, CVTMLE = CVTMLE, byfold_Q = byfold_Q,
      #                           models = models_control, max_iter = max_iter, adapt_stop = adapt_stop,
      #                           adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps,
      #                           return_fW = return_fW, maxpY = maxpY, TMLE_updater = TMLE_updater, verbose = verbose)
      #   return(res)
    #   }
    #   res_byt[] <- res_byt[rev(seq_along(tvals))] # re-assign to order results by increasing t
    # } else {
    #   res_byt <- vector(mode = "list", length = length(tvals))
    #   for (t_idx in rev(seq_along(tvals))) {
    #     t_period <- tvals[t_idx]
    #     res <- fit_GCOMP_onet(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule,
    #                             TMLE = TMLE, iterTMLE = iterTMLE, CVTMLE = CVTMLE, byfold_Q = byfold_Q,
    #                             models = models_control, max_iter = max_iter, adapt_stop = adapt_stop,
    #                             adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps,
    #                             return_fW = return_fW, maxpY = maxpY, TMLE_updater = TMLE_updater, verbose = verbose)
    #     res_byt[[t_idx]] <- res
    #   }
    # }
  # )

  ## Restore backed up nodes, even in the event of failure (otherwise the input data is corrupted for good)
  OData$restoreNodes(c(intervened_TRT,intervened_MONITOR))
  ## remove the newly added rule name column from observed data
  if (added_rule_col) {
    OData$dat.sVar[, "rule.name" := NULL]
  }

  if (inherits(tmle.run.res, "try-error")) {
    stop(
"...attempt at running TMLE for one or several time-points has failed for unknown reason;
If this error cannot be fixed, consider creating a replicable example and filing a bug report at:
  https://github.com/osofr/stremr/issues
", call. = TRUE)
  }

  resultDT <- rbindlist(res_byt)
  ## to extract the EIC estimates by time point (no longer returning those separately, only as part of the estimates data.table)
  # ICs_byt <- resultDT[["IC.St"]]
  # IC.Var.S.d <- t(do.call("cbind", ICs_byt))

  resultDT <- cbind(est_name = est_name, resultDT)
  if (!"rule.name" %in% names(resultDT)) {
    resultDT[, "rule.name" := eval(as.character(rule_name))]
  }

  attr(resultDT, "estimator_short") <- est_name
  attr(resultDT, "estimator_long") <- est_name
  attr(resultDT, "nID") <- OData$nuniqueIDs
  attr(resultDT, "rule_name") <- rule_name
  attr(resultDT, "stratifyQ_by_rule") <- stratifyQ_by_rule
  attr(resultDT, "stratify_by_last") <- stratify_by_last
  attr(resultDT, "trunc_weights") <- trunc_weights
  attr(resultDT, "time") <- resultDT[["time"]]

  res_out <- list(
              # est_name = est_name,
              # periods = tvals,
              # IC.Var.S.d = list(IC.S = IC.Var.S.d),
              # nID = OData$nuniqueIDs,
              wts_data = { if ((TMLE || iterTMLE) && return_wts) { IPWeights } else { NULL } },
              # rule_name = rule_name,
              # trunc_weights = trunc_weights
              estimates = resultDT
              )

  attr(res_out, "estimator_short") <- est_name
  attr(res_out, "estimator_long") <- est_name
  return(res_out)
}

# ------------------------------------------------------------------------------------------------
# ITERATIVE (UNIVARIATE) TMLE AGLORITHM for a single time-point.
# Called as part of the fit_GCOMP_onet()
# ------------------------------------------------------------------------------------------------
iterTMLE_onet <- function(OData, Qlearn.fit, Qreg_idx, max_iter = 15, adapt_stop = TRUE, adapt_stop_factor = 10, tol_eps = 0.001) {
  get_field_Qclass <- function(allQmodels, fieldName) {
    lapply(allQmodels, function(Qclass) Qclass[[fieldName]])
    # lapply(allQmodels, function(Qclass) Qclass$getPsAsW.models()[[1]][[fieldName]])
  }

  eval_idx_used_to_fit_initQ <- function(allQmodels, OData) {
    lapply(allQmodels, function(Qclass) which(Qclass$define_idx_to_fit_initQ(data = OData)))
    # lapply(allQmodels, function(Qclass) which(Qclass$getPsAsW.models()[[1]]$define_idx_to_fit_initQ(data = OData)))
  }

  # clean up left-over indices after we are done iterating
  eval_wipe.all.indices <- function(allQmodels, OData) {
    lapply(allQmodels, function(Qclass) Qclass$wipe.all.indices)
    # lapply(allQmodels, function(Qclass) Qclass$getPsAsW.models()[[1]]$wipe.all.indices)
  }

  Propagate_TMLE_fits <- function(allQmodels, OData, TMLE.fit) {
    lapply(allQmodels,
      function(Qclass) Qclass$Propagate_TMLE_fit(data = OData, new.TMLE.fit = TMLE.fit))
      # function(Qclass) Qclass$getPsAsW.models()[[1]]$Propagate_TMLE_fit(data = OData, new.TMLE.fit = TMLE.fit))
    return(invisible(allQmodels))
  }

  allQmodels <- Qlearn.fit$getPsAsW.models() # Get the individual Qlearning classes
  # res_all_subset_idx <- as.vector(sort(unlist(get_field_Qclass(allQmodels, "subset_idx"))))
  # use_subset_idx <- res_all_subset_idx
  # use_subset_idx <- res_idx_used_to_fit_initQ

  # one cat'ed vector of all observations that were used for fitting init Q & updating TMLE (across all t's):
  res_idx_used_to_fit_initQ <- as.vector(sort(unlist(get_field_Qclass(allQmodels, "idx_used_to_fit_initQ"))))
  # res_idx_used_to_fit_initQ_2 <- as.vector(sort(unlist(eval_idx_used_to_fit_initQ(allQmodels, OData))))
  # all.equal(res_idx_used_to_fit_initQ, res_idx_used_to_fit_initQ_2)

  # Consider only observations with non-zero weights, these are the only obs that are needed for the TMLE update:
  idx_all_wts_above0 <- which(OData$IPwts_by_regimen[["cum.IPAW"]] > 0)
  use_subset_idx <- intersect(idx_all_wts_above0, res_idx_used_to_fit_initQ)
  wts_TMLE <- OData$IPwts_by_regimen[use_subset_idx, "cum.IPAW", with = FALSE][[1]]

  for (iter in 1:(max_iter+1)) {
    Qkplus1 <- OData$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]
    Qk_hat <- OData$dat.sVar[use_subset_idx, "Qk_hat", with = FALSE][[1]]

    # ------------------------------------------------------------------------------------
    # ESTIMATE OF THE EIC:
    # ------------------------------------------------------------------------------------
    # Get t-specific and i-specific components of the EIC for all t > t.init:
    EIC_i_tplus <- wts_TMLE * (Qkplus1 - Qk_hat)
    # Get t-specific and i-specific components of the EIC for all t = t.init (mean pred from last reg, all n obs)
    res_lastPredQ <- allQmodels[[Qreg_idx[1]]]$predictAeqa()  # Qreg_idx[1] is the index for the last Q-fit
    # res_lastPredQ <- allQmodels[[length(allQmodels)]]$predictAeqa()

    EIC_i_t0 <- (res_lastPredQ - mean(res_lastPredQ))
    # Sum them all up and divide by N -> obtain the estimate of P_n(D^*_n):
    EIC_est <- (sum(EIC_i_t0) + sum(EIC_i_tplus)) / OData$nuniqueIDs
    if (gvars$verbose) print("Mean estimate of the EIC: " %+% EIC_est)

    # Quit loop if the error tolerance level has been reached
    if (iter > 1) {
      if (adapt_stop) {
        # if (abs(EIC_est) < (1 / OData$nuniqueIDs)) break
        if (abs(EIC_est) < (1 / (adapt_stop_factor*sqrt(OData$nuniqueIDs)))) break
      } else {
        if (!is.null(tol_eps) & (abs(TMLE.fit$TMLE_intercept) <= tol_eps)) break
      }
    } else if (iter > max_iter) {
      break
    }

    TMLE.fit <- tmle.update(Qkplus1 = Qkplus1, Qk_hat = Qk_hat, IPWts = wts_TMLE, lower_bound_zero_Q = FALSE, skip_update_zero_Q = FALSE)
    Propagate_TMLE_fits(allQmodels, OData, TMLE.fit)
  }

  if (gvars$verbose) print("iterative TMLE ran for N iter: " %+% iter)

  # clean up after we are done iterating (set the indices in Q classes to NULL to conserve memory)
  tmp <- eval_wipe.all.indices(allQmodels)

  # EVALUTE THE t and i-specific components of the EIC (estimates):
  # Qkplus1 <- OData$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]
  # Qk_hat <- OData$dat.sVar[use_subset_idx, "Qkplus1", with = FALSE][[1]]
  # EIC_i_t_calc <- wts_TMLE * (Qkplus1 - Qk_hat)
  # OData$dat.sVar[use_subset_idx, ("EIC_i_t") := EIC_i_t_calc]
  OData$dat.sVar[use_subset_idx, ("EIC_i_t") := EIC_i_tplus]

  return(invisible(Qlearn.fit))
}



