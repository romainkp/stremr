# ------------------------------------------------------------------------------------------
# H2O ISSUES WITH PARALLEL TRAINING:
# ------------------------------------------------------------------------------------------
# 1) EVERYTIME a subset [,] is done on frame (H2O.dat.sVar[rows_subset, vars]) a new frame is created with an automatically ID
# 2) EVERYTIME a prediction is made a new h2o frame is created: predictions_"modelIDname"_on_"H2OFRAMENAME"
# 3) EVERYTIME a prediction result is pulled ([,"p1"]), a new temp FRAME  is created
# 4) WHEN DOING "subsetH2Oframe[, outvar]" in BinomialH2O$setdata A NEW temp FRAME (AUTOMATIC ID) is CREATED
# 5) Having Q.kplus1 as a column in the main DT or FRAME is also an issue, since parallel training might overwrite it

# ------------------------------------------------------------------------------------------
# POTENTIAL SOLUTION:
# ------------------------------------------------------------------------------------------
# *) PASS Q.kplus1 AS a vector argument?
# *) using existing approach to LOADING A NEW FRAME FOR EVERY single call to QlearnModel$fit,
#    but pass the args so that destination_frame is always different for each separate call of fitSeqGcomp_onet()
# *) remove all the implicit FRAME creating steps ('[,]' with automatic IDs) with explicit steps so that destination_frame can be specified
# *) use h2o.removeAll() to remove all frames at once
# *) use h2o.rm(ids) to remove a single frame by its ID
# *) h2o.ls() to list all current frames
# *) Classify a single instance at a time:
# http://www.h2o.ai/product/faq/#H2OClassifyInstance
# The plain Java (POJO) scoring predict API is: public final float[] predict( double[] data, float[] preds)
# // Pass in data in a double[], pre-aligned to the Model's requirements.
# // Jam predictions into the preds[] array; preds[0] is reserved for the
# // main prediction (class for classifiers or value for regression),
# // and remaining columns hold a probability distribution for classifiers.

# ------------------------------------------------------------------------------------------
# TO DO:
# ------------------------------------------------------------------------------------------
# *) Need to finish/clean-up the suit of tests for all componensts.
# *) BUG with TMLE weights when applying to real data -> last row with NA's screws everything up.
# *) The definition of Qperiods below needs to be based on actual periods observed in the data.
# *) Consider bringing in tmlenet syntax for defining the interventions in a node-like style.
# *) Need to be able to differentiate binomial (binary) outcome classification and continuous outcome regression for h2o:
#    1. Setting / not setting the outcome as factor
#    2. Setting the GBM distribution to "bernoulli"/"gaussian"
#    3. Validating the h2o.glm with bernoulli and continous outcome, setting the solver to IRLSM
# ------------------------------------------------------------------------------------------
# (DONE) Need to make sure all functions in main_estimation.R work with the modified class structure.
# (DONE) Redo getIPWeights() to work with counterfactual A's and N's rather than rule followers.
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
# Outstanding issues:
# ------------------------------------------------------------------------------------------
# *** Stochastic interventions ***
# *) intervened_TRT/intervened_MONITOR are vectors of counterfactual probs -> need to allow each to be multivariate (>1 cols)
# *) For column 0<intervened_TRT<1 it defines the counterfactual probability that P(TRT[t]^*=1)=intervened_TRT[t].
# *) gstar can be a vector, i.e., abar=(0,0,0,0) or a matrix of counterfactual treatments, like in ltmle.
#     Deal with stochastic interventions in QlearnModel class (possibly daughter classes) with direct intergration
#     eval'ed as direct weighted sum of predicted Q's with weights given by g.star (on A and N).
# ------------------------------------------------------------------------------------------
#  *** rule_followers:
#     Once you go off the treatment first time, this is it, the person is censored for the rest of the follow-up.
#     Rule followers are now evaluated automatically by comparing (intervened_TRT and TRT and CENS) and (intervened_MONITOR and MONITOR and CENS).
#     Rule followers are: (intervened_TRT = 1) & (TRT == 1) or (intervened_TRT = 0) & (TRT == 0) or (intervened_TRT > 0 & intervened_TRT < 0) & (Not Censored).
#     Exactly the same logic is also applied to (intervened_MONITOR & MONITOR) when these are specified.
#     NOT TESTED: If either TRT or MONITOR is multivariate (>1 col), this logic needs to be applied FOR EACH COLUMN of TRT/MONITOR.
# ------------------------------------------------------------------------------------------
# *** Accessing correct QlearnModel ***
#     Need to come up with a good way of accessing correct terminal QlearnModel to get final Q-predictions in private$probAeqa.
#     These predictions are for the final n observations that were used for evaluating E[Y^a].
#     However, if we do stratify=TRUE and use a stack version of g-comp with different strata all in one stacked database it will be difficult
#     to pull the right QlearnModel.
#     Alaternative is to ignore that completely and go directly for the observed data (by looking up the right column/row combos).
# ------------------------------------------------------------------------------------------
# *** Pooling & performing only one TMLE update across all Q.k ***
#     In general, we can pool all Q.kplus (across all t's) and do a single update on all of them
#     Would make this TMLE iterative, but would not require refitting on the initial Q's
#     See how fast is the G-comp formula prediction without fitting (the subseting/prediction loop alone) -> maybe fast enough to allow iterations
# ------------------------------------------------------------------------------------------
# *** Pooling across regimens ***
#     Since new regimen results in new Q.kplus1 and hence new outcome -> requires a separate regression for each regimen:
#     => Can either pool all Q.regimens at once (fit one model for repated observations, one for each rule, smoothing over rules).
#     => Can use the same stacked dataset of regime-followers, but with a separate stratification for each regimen.
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
# Sequential (Recursive) G-COMP Algorithm:
# ------------------------------------------------------------------------------------------
# 1. At each t iteration:
#     Fitting:
#        Outcome is always Q.kplus1 (either by regimen or not). Fitting subset excludes all obs censored at t.
#     Prediction:
#        Swap the entire column(s) A[t] with A^*[t] by renaming them in OData$dat.sVar (might do this for more than one regimen in the future, if pooling Q's)
#        Subset all observation by t (including all censored and non-followers).
#        PredictP1() for everybody in the subset, then save the prediction in rows Q.kplus1[t] and rows Q.kplus1[t-1] for all obs that were used in prediction.
#        Results in correct algorithm even when stratifyQ_by_rule=TRUE: next iteration (t-1) will do fit only based on obs that followed the rule at t-1.
# 2. At next iteration t-1:
#     Fitting:
#        Use the predictions in Q.kplus1[t-1] as new outcomes and repeat until reached minimum t.
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# When useonly_t_TRT or useonly_t_MONITOR is specified, need to set nodes to their observed values, rather than the counterfactual values
# ------------------------------------------------------------------------------------------------------------------------
# Do it separately for gstar_TRT & gstar_MONITOR
# Loop over each node in gstar_TRT / gstar_MONITOR
# Do it only once for all observations inside main tmle call
# Back-up a copy of all gstar nodes first, the original copy is then restored when finished running
# The observations which get swapped with g0 values are defined by:
# subset_idx <- OData$evalsubst(subset_exprs = useonly_t_NODE)
# probability of P(A^*(t)=n(t)) or P(N^*(t)=n(t)) under counterfactual A^*(t) or N^*(t) and observed a(t) or n(t)
# Example call:
# defineNodeGstarGComp(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, stratifyQ_by_rule)
defineNodeGstarGComp <- function(OData, intervened_NODE, NodeNames, useonly_t_NODE, stratifyQ_by_rule) {
  # if intervened_NODE returns more than one rule-column, evaluate g^* for each and the multiply to get a single joint (for each time point)
  if (!is.null(intervened_NODE)) {
    gstar.NODEs <- intervened_NODE
    for (intervened_NODE_col in intervened_NODE) CheckVarNameExists(OData$dat.sVar, intervened_NODE_col)
    assert_that(length(intervened_NODE) == length(NodeNames))

    # ------------------------------------------------------------------------------------------
    # Modify the observed input intervened_NODE in OData$dat.sVar with values from NodeNames for subset_idx:
    # ------------------------------------------------------------------------------------------
    subset_idx <- OData$evalsubst(subset_exprs = useonly_t_NODE)
    OData$replaceNodesVals(!subset_idx, nodes_to_repl = intervened_NODE, source_for_repl = NodeNames)
    # ------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------
    # update rule followers for trt if doing stratified G-COMP:
    # Note this will define rule followers based on REPLACED intervened_NODE in dat.sVar (i.e., modified n^*(t) under N.D.E.)
    # FOR NDE BASED TMLE THE DEFINITION OF RULE-FOLLOWERS CHANGES ACCORDINGLY based on modified n^*(t) and a^*(t)
    # ------------------------------------------------------------------------------------------
    if (stratifyQ_by_rule) {
      rule_followers_idx <- OData$eval_rule_followers(NodeName = NodeNames, gstar.NodeName = intervened_NODE)
      OData$rule_followers_idx <- rule_followers_idx & OData$rule_followers_idx & OData$uncensored_idx
    }

  } else {
    # use the actual (observed) node names under g0:
    gstar.NODEs <- NodeNames
  }
  return(gstar.NODEs)
}

# ---------------------------------------------------------------------------------------
#' Iterative TMLE wrapper for \code{fitSeqGcomp}
#'
#' Calls \code{fitSeqGcomp} with argument \code{iterTMLE = TRUE}.
#' @param ... Arguments that will be passed down to the underlying function \code{fitSeqGcomp}
#' @return \code{data.table} with survival by time for sequential GCOMP and iterative TMLE
#' @seealso \code{\link{fitSeqGcomp}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fitIterTMLE <- function(...) {
  fitSeqGcomp(TMLE = FALSE, iterTMLE = TRUE, ...)
}

# ---------------------------------------------------------------------------------------
#' TMLE wrapper for \code{fitSeqGcomp}
#'
#' Calls \code{fitSeqGcomp} with argument \code{TMLE = TRUE}.
#' @param ... Arguments that will be passed down to the underlying function \code{fitSeqGcomp}
#' @return \code{data.table} with TMLE survival by time
#' @seealso \code{\link{fitSeqGcomp}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fitTMLE <- function(...) {
  fitSeqGcomp(TMLE = TRUE, ...)
}

# ---------------------------------------------------------------------------------------
#' Fit sequential GCOMP and TMLE for survival
#'
#' Interventions on up to 3 nodes are allowed: \code{CENS}, \code{TRT} and \code{MONITOR}.
#' TMLE adjustment will be based on the inverse of the propensity score fits for the observed likelihood (g0.C, g0.A, g0.N),
#' multiplied by the indicator of not being censored and the probability of each intervention in \code{intervened_TRT} and \code{intervened_MONITOR}.
#' Requires column name(s) that specify the counterfactual node values or the counterfactual probabilities of each node being 1 (for stochastic interventions).
#' @param OData Input data object created by \code{importData} function.
#' @param t_periods Specify the vector of time-points for which the survival function (and risk) should be estimated
#' @param Qforms Regression formulas, one formula per Q. Only main-terms are allowed.
#' @param Qstratify Placeholder for future user-defined model stratification for all Qs (CURRENTLY NOT FUNCTIONAL, WILL RESULT IN ERROR)
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
#' @param stratifyQ_by_rule Set to \code{TRUE} for stratification, fits the outcome model (Q-learning) among rule-followers only.
#' Setting to \code{FALSE} will fit the outcome model (Q-learning) across all observations (pooled regression).
#' @param TMLE Set to \code{TRUE} to run the usual longitudinal TMLE algorithm (with a separate TMLE update of Q for every sequential regression).
#' @param iterTMLE Set to \code{TRUE} to run the iterative univariate TMLE instead of the usual longitudinal TMLE.
#' When set to \code{TRUE} this will also provide the standard sequential Gcomp as party of the output.
#' Must set \code{TMLE}=\code{FALSE} when setting this to \code{TRUE}.
#' @param IPWeights (Optional) result of calling function \code{getIPWeights} for running TMLE (evaluated automatically when missing)
#' @param stabilize Set to \code{TRUE} to use stabilized weights for the TMLE
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @param params_Q Optional parameters to be passed to the specific fitting algorithm for Q-learning
#' @param weights Optional \code{data.table} with additional observation-time-specific weights.  Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param max_iter For iterative TMLE only: Integer, set to maximum number of iterations for iterative TMLE algorithm.
#' @param adapt_stop For iterative TMLE only: Choose between two stopping criteria for iterative TMLE, default is \code{TRUE},
#' which will stop the iterative TMLE algorithm in an adaptive way. Specifically, the iterations will stop when the mean estimate
#' of the efficient influence curve is less than or equal to 1 / (\code{adapt_stop_factor}*sqrt(\code{N})), where
#' N is the total number of unique subjects in data and \code{adapt_stop_factor} is set to 10 by default.
#' When \code{TRUE}, the argument \code{tol_eps} is ignored and TMLE stops when either \code{max_iter} has been reached or this criteria has been satisfied.
#' When \code{FALSE}, the stopping criteria is determined by values of \code{max_iter} and \code{tol_eps}.
#' @param adapt_stop_factor For iterative TMLE only: The adaptive factor to choose the stopping criteria for iterative TMLE when \code{adapt_stop} is set to \code{TRUE}. Default is 10.
#' TMLE will keep iterative until
#' the mean estimate of the efficient influence curve is less than 1 / (\code{adapt_stop_factor}*sqrt(\code{N})) or when the number of iterations is \code{max_iter}.
#' @param tol_eps For iterative TMLE only: Numeric error tolerance for the iterative TMLE update.
#' The iterative TMLE algorithm will stop when the absolute value of the TMLE intercept update is below \code{tol_eps}
#' @param parallel Set to \code{TRUE} to run the sequential Gcomp or TMLE in parallel (uses \code{foreach} with \code{dopar} and requires a previously defined parallel back-end cluster)
#' @param verbose ...
#' @return ...
#' @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
fitSeqGcomp <- function(OData, t_periods,
                        Qforms, Qstratify = NULL,
                        intervened_TRT = NULL, intervened_MONITOR = NULL,
                        useonly_t_TRT = NULL, useonly_t_MONITOR = NULL,
                        rule_name = paste0(c(intervened_TRT, intervened_MONITOR), collapse = ""),
                        stratifyQ_by_rule = FALSE,
                        TMLE = FALSE,
                        iterTMLE = FALSE,
                        IPWeights = NULL,
                        stabilize = FALSE,
                        trunc_weights = 10^6,
                        params_Q = list(),
                        weights = NULL,
                        max_iter = 15,
                        adapt_stop = TRUE,
                        adapt_stop_factor = 10,
                        tol_eps = 0.001,
                        parallel = FALSE,
                        verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  assert_that(is.list(params_Q))
  assert_that(is.logical(adapt_stop))

  if (TMLE & iterTMLE) stop("Either 'TMLE' or 'iterTMLE' must be set to FALSE. Cannot estimate both within a single algorithm run.")

  if (missing(rule_name)) rule_name <- paste0(c(intervened_TRT,intervened_MONITOR), collapse = "")
  # ------------------------------------------------------------------------------------------------
  # **** Evaluate the uncensored and initialize rule followers (everybody is a follower by default)
  # ------------------------------------------------------------------------------------------------

  OData$uncensored_idx <- OData$eval_uncensored()
  OData$rule_followers_idx <- rep.int(TRUE, nrow(OData$dat.sVar)) # (everybody is a follower by default)

  # ------------------------------------------------------------------------------------------------
  # **** Load dat.sVar into H2O.Frame memory if loading data only once
  # ------------------------------------------------------------------------------------------------
  # mainH2Oframe <- OData$fast.load.to.H2O(OData$dat.sVar,
  #                                       saveH2O = TRUE,
  #                                       destination_frame = "H2OMainDataTable")
  # ------------------------------------------------------------------------------------------------
  # *** TO change the column names of H2O.FRAME (for Qlearning steps if data is loaded only once) ***
  # ------------------------------------------------------------------------------------------------
  # h2o::colnames(subsetH2Oframe) <- c("T", "C", "h", "N")
  # names(subsetH2Oframe) <- names(subsetH2Oframe)%+%"_1"
  # names(subsetH2Oframe)[1] <- "T2"%+%"_1"

  # ------------------------------------------------------------------------------------------------
  # **** Define regression paramers specific to Q-learning (continuous outcome)
  # ------------------------------------------------------------------------------------------------
  # parameter for running h2o.glm with continous outcome (this is the only sovler that works)
  # to try experimental solvers in h2o.glm (see line #901 of GLM.java)
  # params_Q$solver <- "COORDINATE_DESCENT"
  # params_Q$solver <- "COORDINATE_DESCENT_NAIVE"
  # params_Q$solver <- "IRLSM"
  # params_Q$solver = "L_BFGS"

  # parameter for running GBM with continuous outcome (default is classification with "bernoulli" and 0/1 outcome):
  params_Q$distribution <- "gaussian"

  # ------------------------------------------------------------------------------------------------
  # **** Add weights if TMLE=TRUE and if weights were defined
  # NOTE: This needs to be done only once if evaluating survival over several t_periods
  # ------------------------------------------------------------------------------------------------
  if (TMLE || iterTMLE) {
    if (is.null(IPWeights)) {
      if (gvars$verbose) message("...evaluating IPWeights for TMLE...")
      IPWeights <- getIPWeights(OData, intervened_TRT, intervened_MONITOR, useonly_t_TRT, useonly_t_MONITOR, rule_name)
      if (stabilize) IPWeights[, "cum.IPAW" := eval(as.name("cum.stab.P")) * eval(as.name("cum.IPAW"))]
      if (trunc_weights < Inf) IPWeights[eval(as.name("cum.IPAW")) > trunc_weights, ("cum.IPAW") := trunc_weights]
    } else {
      getIPWeights_fun_call <- attributes(IPWeights)[['getIPWeights_fun_call']]
      if (gvars$verbose) message("applying user-specified IPWeights obtained with a call: \n" %+% deparse(getIPWeights_fun_call)[[1]])
      assert_that(all.equal(attributes(IPWeights)[['intervened_TRT']], intervened_TRT))
      assert_that(all.equal(attributes(IPWeights)[['intervened_MONITOR']], intervened_MONITOR))
      # assert_that(attributes(IPWeights)[['stabilize']] == FALSE)
    }
    assert_that(is.data.table(IPWeights))
    assert_that("cum.IPAW" %in% names(IPWeights))
    OData$IPwts_by_regimen <- IPWeights

    if (!is.null(weights)) stop("optional argument 'weights' is not implemented for TMLE or GCOMP")
    # Add additional observation-specific weights to the cumulative weights:
    # IPWeights <- process_opt_wts(IPWeights, weights, nodes, adjust_outcome = FALSE)
  }

  if (missing(t_periods)) stop("must specify survival 't_periods' of interest (time period values from column " %+% nodes$tnode %+% ")")

  # ------------------------------------------------------------------------------------------
  # Create a back-up of the observed input gstar nodes (created by user in input data):
  # Will add new columns (that were not backed up yet) TO SAME backup data.table
  # ------------------------------------------------------------------------------------------
  OData$backupNodes(c(intervened_TRT,intervened_MONITOR))

  # ------------------------------------------------------------------------------------------------
  # Define the intervention nodes
  # Modify the observed input intervened_NODE in OData$dat.sVar with values from NodeNames for subset_idx
  # ------------------------------------------------------------------------------------------------
  # browser()

  gstar.A <- defineNodeGstarGComp(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, stratifyQ_by_rule)
  gstar.N <- defineNodeGstarGComp(OData, intervened_MONITOR, nodes$Nnodes, useonly_t_MONITOR, stratifyQ_by_rule)
  interventionNodes.g0 <- c(nodes$Anodes, nodes$Nnodes)
  interventionNodes.gstar <- c(gstar.A, gstar.N)

  # When the same node names belongs to both g0 and gstar it doesn't need to be intervened upon, so exclude
  common_names <- intersect(interventionNodes.g0, interventionNodes.gstar)
  interventionNodes.g0 <- interventionNodes.g0[!interventionNodes.g0 %in% common_names]
  interventionNodes.gstar <- interventionNodes.gstar[!interventionNodes.gstar %in% common_names]

  OData$interventionNodes.g0 <- interventionNodes.g0
  OData$interventionNodes.gstar <- interventionNodes.gstar

  # ------------------------------------------------------------------------------------------------
  # RUN GCOMP OR TMLE FOR SEVERAL TIME-POINTS EITHER IN PARALLEL OR SEQUENTIALLY
  # ------------------------------------------------------------------------------------------------
  est_name <- ifelse(TMLE, "TMLE", ifelse(iterTMLE, "GCOMP & Iter.TMLE", "GCOMP"))
  tmle.run.res <- try(
    if (parallel) {
      mcoptions <- list(preschedule = FALSE)
      '%dopar%' <- foreach::'%dopar%'
      res_byt <- foreach::foreach(t_idx = seq_along(t_periods), .options.multicore = mcoptions) %dopar% {
        t_period <- t_periods[t_idx]
        res <- fitSeqGcomp_onet(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule, TMLE = TMLE, iterTMLE = iterTMLE,
                                params_Q = params_Q, max_iter = max_iter, adapt_stop = adapt_stop, adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps, verbose = verbose)
        return(res)
      }
    } else {
      res_byt <- vector(mode = "list", length = length(t_periods))
      for (t_idx in seq_along(t_periods)) {
        t_period <- t_periods[t_idx]
        res <- fitSeqGcomp_onet(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule, TMLE = TMLE, iterTMLE = iterTMLE,
                                params_Q = params_Q, max_iter = max_iter, adapt_stop = adapt_stop, adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps, verbose = verbose)
        res_byt[[t_idx]] <- res
      }
    }
  )

  # ------------------------------------------------------------------------------------------------
  # Restore backed up nodes, even in the event of failure (otherwise the input data is corrupted for good)
  # ------------------------------------------------------------------------------------------------
  OData$restoreNodes(c(intervened_TRT,intervened_MONITOR))

  if (inherits(tmle.run.res, "try-error")) { # TMLE update failed
    stop(
"...attempt at running TMLE for one or several time-points has failed for unknown reason;
If this error cannot be fixed, consider creating a replicable example and filing a bug report at:
  https://github.com/osofr/stremr/issues
", call. = TRUE)
  }

  ICs_byt <- lapply(res_byt, '[[', "IC_i_onet")
  IC.Var.S.d <- t(do.call("cbind", ICs_byt))
  # IC.Var.S.d <- matrix(NA, nrow = length(ICs_byt), ncol = length(ICs_byt[[1]]))

  ests_byt <- lapply(res_byt, '[[', "resDF_onet")
  resultDT <- data.table(est_name = est_name, rbindlist(ests_byt))
  resultDT[, "rule.name" := eval(as.character(rule_name))]

  res_out <- list(
              estimates = resultDT,
              est_name = est_name,
              periods = t_periods,
              IC.Var.S.d = list(IC.S = IC.Var.S.d),
              nID = OData$nuniqueIDs,
              wts_data = { if (TMLE || iterTMLE) {IPWeights} else {NULL}},
              rule_name = rule_name,
              trunc_weights = trunc_weights)

  return(res_out)
}

# ------------------------------------------------------------------------------------------------
# ITERATIVE (UNIVARIATE) TMLE AGLORITHM for a single time-point.
# Called as part of the fitSeqGcomp_onet()
# ------------------------------------------------------------------------------------------------
iterTMLE_onet <- function(OData, Qlearn.fit, Qreg_idx, max_iter = 15, adapt_stop = TRUE, adapt_stop_factor = 10, tol_eps = 0.001) {
  get_field_Qclass <- function(allQmodels, fieldName) {
    lapply(allQmodels, function(Qclass) Qclass$getPsAsW.models()[[1]][[fieldName]])
  }

  eval_idx_used_to_fit_initQ <- function(allQmodels, OData) {
    lapply(allQmodels, function(Qclass) which(Qclass$getPsAsW.models()[[1]]$define_idx_to_fit_initQ(data = OData)))
  }

  # clean up left-over indices after we are done iterating
  eval_wipe.all.indices <- function(allQmodels, OData) {
    lapply(allQmodels, function(Qclass) Qclass$getPsAsW.models()[[1]]$wipe.all.indices)
  }

  Propagate_TMLE_fits <- function(allQmodels, OData, TMLE.fit) {
    lapply(allQmodels,
      function(Qclass) Qclass$getPsAsW.models()[[1]]$Propagate_TMLE_fit(data = OData, new.TMLE.fit = TMLE.fit))
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
    prev_Q.kplus1 <- OData$dat.sVar[use_subset_idx, "prev_Q.kplus1", with = FALSE][[1]]
    init_Q_fitted_only <- OData$dat.sVar[use_subset_idx, "Q.kplus1", with = FALSE][[1]]

    # ------------------------------------------------------------------------------------
    # ESTIMATE OF THE EIC:
    # ------------------------------------------------------------------------------------
    # Get t-specific and i-specific components of the EIC for all t > t.init:
    EIC_i_tplus <- wts_TMLE * (prev_Q.kplus1 - init_Q_fitted_only)
    # Get t-specific and i-specific components of the EIC for all t = t.init (mean pred from last reg, all n obs)
    res_lastPredQ <- Qlearn.fit$predictRegK(Qreg_idx[1], OData$nuniqueIDs) # Qreg_idx[1] is the index for the last Q-fit
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
        if (!is.null(tol_eps) & (abs(TMLE.fit$TMLE.intercept) <= tol_eps)) break
      }
    } else if (iter > max_iter) {
      break
    }

    TMLE.fit <- tmle.update(prev_Q.kplus1 = prev_Q.kplus1, init_Q_fitted_only = init_Q_fitted_only, IPWts = wts_TMLE, lower_bound_zero_Q = FALSE, skip_update_zero_Q = FALSE)
    Propagate_TMLE_fits(allQmodels, OData, TMLE.fit)
  }

  if (gvars$verbose) print("iterative TMLE ran for N iter: " %+% iter)

  # clean up after we are done iterating (set the indices in Q classes to NULL to conserve memory)
  tmp <- eval_wipe.all.indices(allQmodels)

  # EVALUTE THE t-specific and i-specific components of the EIC (estimates):
  # prev_Q.kplus1 <- OData$dat.sVar[use_subset_idx, "prev_Q.kplus1", with = FALSE][[1]]
  # init_Q_fitted_only <- OData$dat.sVar[use_subset_idx, "Q.kplus1", with = FALSE][[1]]
  # EIC_i_t_calc <- wts_TMLE * (prev_Q.kplus1 - init_Q_fitted_only)
  # OData$dat.sVar[use_subset_idx, ("EIC_i_t") := EIC_i_t_calc]
  OData$dat.sVar[use_subset_idx, ("EIC_i_t") := EIC_i_tplus]

  return(invisible(Qlearn.fit))
}

fitSeqGcomp_onet <- function(OData, t_period, Qforms, Qstratify, stratifyQ_by_rule, TMLE, iterTMLE, params_Q, max_iter = 10, adapt_stop = TRUE, adapt_stop_factor = 10, tol_eps = 0.001,
                             verbose = getOption("stremr.verbose")) {
  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names

  # browser()
  # ------------------------------------------------------------------------------------------------
  # Defining the t periods to loop over FOR A SINGLE RUN OF THE iterative G-COMP/TMLE (one survival point)
  # **** TO DO: The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------
  Qperiods <- rev(OData$min.t:t_period)
  Qreg_idx <- rev(seq_along(Qperiods))
  Qstratas_by_t <- as.list(nodes[['tnode']] %+% " == " %+% (Qperiods))
  names(Qstratas_by_t) <- rep.int("Q.kplus1", length(Qstratas_by_t))
  # Adding user-specified stratas to each t Q-regression:
  all_Q_stratify <- Qstratas_by_t
  if (!is.null(Qstratify)) {
    assert_that(is.vector(Qstratify))
    assert_that(is.character(Qstratify))
    for (idx in seq_along(all_Q_stratify)) {
      all_Q_stratify[[idx]] <- stringr::str_c(all_Q_stratify[[idx]], " & ", Qstratify)
    }
  }

  # ------------------------------------------------------------------------------------------------
  # **** Process the input formulas and stratification settings
  # **** TO DO: Add checks that Qforms has correct number of regressions in it
  # ------------------------------------------------------------------------------------------------
  Qforms.default <- rep.int("Q.kplus1 ~ Lnodes + Anodes + Cnodes + Nnodes", length(Qperiods))
  if (missing(Qforms)) {
    Qforms_single_t <- Qforms.default
  } else {
    Qforms_single_t <- Qforms[seq_along(Qperiods)]
  }

  # ------------------------------------------------------------------------------------------------
  #  ****** Q.kplus1 THIS NEEDS TO BE MOVED OUT OF THE MAIN data.table *****
  # Running fitSeqGcomp_onet might conflict with different Q.kplus1
  # ------------------------------------------------------------------------------------------------
  # **** G-COMP: Initiate Q.kplus1 - (could be multiple if more than one regimen)
  # That column keeps the tabs on the running Q-fit (SEQ G-COMP)
  # ------------------------------------------------------------------------------------------------
  # set the initial (default values of the t-specific and i-specific EIC estimates):
  OData$dat.sVar[, ("EIC_i_t") := 0.0]
  # set the initial values of Q (the observed outcome node):
  OData$dat.sVar[, "Q.kplus1" := as.numeric(get(OData$nodes$Ynode))]

  OData$set.sVar.type(name.sVar = "Q.kplus1", new.type = "binary")
  OData$set.sVar.type(name.sVar = "EIC_i_t", new.type = "binary")
  # OData$def.types.sVar() # bottleneck

  # ------------------------------------------------------------------------------------------------
  # **** Define regression classes for Q.Y and put them in a single list of regressions.
  # **** TO DO: This could also be done only once in the main routine, then just subset the appropriate Q_regs_list
  # ------------------------------------------------------------------------------------------------
  Q_regs_list <- vector(mode = "list", length = length(Qstratas_by_t))
  names(Q_regs_list) <- unlist(Qstratas_by_t)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")
  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms_single_t[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    reg <- RegressionClassQlearn$new(Qreg_counter = Qreg_idx[i], t_period = Qperiods[i],
                                     TMLE = TMLE, stratifyQ_by_rule = stratifyQ_by_rule,
                                     outvar = "Q.kplus1", predvars = regform$predvars, outvar.class = list("Qlearn"),
                                     subset_vars = list("Q.kplus1"), subset_exprs = all_Q_stratify[i], model_contrl = params_Q,
                                     censoring = FALSE)
    Q_regs_list[[i]] <- reg
  }

  Qlearn.fit <- GenericModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)

  # Run all Q-learning regressions (one for each subsets defined above, predictions of the last regression form the outcomes for the next:
  Qlearn.fit$fit(data = OData, iterTMLE = iterTMLE)
  OData$Qlearn.fit <- Qlearn.fit

  # When the model is fit with user-defined stratas, need special functions to extract the final fit (current aproach will not work):
  # allQmodels[[1]]$getPsAsW.models()[[1]]$getPsAsW.models()

  # get the individual TMLE updates and evaluate if any updates have failed
  allQmodels <- Qlearn.fit$getPsAsW.models()
  allTMLEfits <- lapply(allQmodels, function(Qmod) Qmod$getPsAsW.models()[[1]]$getTMLEfit)

  TMLEfits <- unlist(allTMLEfits)
  successTMLEupdates <- !is.na(TMLEfits) & !is.nan(TMLEfits)
  ALLsuccessTMLE <- all(successTMLEupdates)
  nFailedUpdates <- sum(!successTMLEupdates)


  # 1a. Grab last reg predictions from Q-regression objects:
  lastQ_inx <- Qreg_idx[1] # The index for the last Q-fit (first time-point)
  # Get the previously saved mean prediction for Q from the very last regression (first time-point, all n obs):
  res_lastPredQ <- Qlearn.fit$predictRegK(lastQ_inx, OData$nuniqueIDs)
  mean_est_t <- mean(res_lastPredQ)
  if (gvars$verbose) print("Surv est: " %+% (1-mean_est_t))
  # # 1b. Can instead grab it directly from the data, using the appropriate strata-subsetting expression
  #   subset_vars <- lastQ.fit$subset_vars
  #   subset_exprs <- lastQ.fit$subset_exprs
  #   subset_idx <- OData$evalsubst(subset_vars = subset_vars, subset_exprs = subset_exprs)
  #   mean(OData$dat.sVar[subset_idx, ][["Q.kplus1"]])

  if (gvars$verbose) print("No. of obs for last prediction of Q: " %+% length(res_lastPredQ))
  if (gvars$verbose) print("EY^* estimate at t="%+%t_period %+%": " %+% round(mean_est_t, 5))

  resDF_onet <- data.frame(t = t_period,
                    risk = mean_est_t,
                    surv = 1 - mean_est_t,
                    ALLsuccessTMLE = ALLsuccessTMLE,
                    nFailedUpdates = nFailedUpdates,
                    type = ifelse(stratifyQ_by_rule, "stratified", "pooled")
                    )

  # ------------------------------------------------------------------------------------------------
  # RUN ITERATIVE TMLE (updating all Q's at once):
  # ------------------------------------------------------------------------------------------------
  if (iterTMLE){
    iter.time <- system.time(
      res <- iterTMLE_onet(OData, Qlearn.fit, Qreg_idx, max_iter = max_iter, adapt_stop = adapt_stop, adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps)
    )
    if (gvars$verbose) {print("Time to run iterative TMLE: "); print(iter.time)}
    # Grab the mean prediction from the very last regression (over all n observations);
    res_lastPredQ <- Qlearn.fit$predictRegK(Qreg_idx[1], OData$nuniqueIDs)
    mean_est_t <- mean(res_lastPredQ)
    if (gvars$verbose) print("Iterative TMLE surv estimate: " %+% (1 - mean_est_t))
    # # # 1b. Grab it directly from the data, using the appropriate strata-subsetting expression
    # lastQ.fit <- Qlearn.fit$getPsAsW.models()[[Qreg_idx[1]]]$getPsAsW.models()[[1]]
    # subset_idx <- OData$evalsubst(subset_vars = lastQ.fit$subset_vars, subset_exprs = lastQ.fit$subset_exprs)
    # res_lastPredQ <- OData$dat.sVar[subset_idx, ][["Q.kplus1"]]
    # mean_est_t  <- mean(res_lastPredQ)
    # print("TMLE surv estimate 2: " %+% (1 - mean_est_t))
    resDF_onet <- cbind(resDF_onet, iterTMLErisk = mean_est_t, iterTMLEsurv = (1 - mean_est_t))
  }

  # ------------------------------------------------------------------------------------------------
  # TMLE INFERENCE
  # ------------------------------------------------------------------------------------------------
  IC_i_onet <- vector(mode = "numeric", length = OData$nuniqueIDs)
  IC_i_onet[] <- NA

  if (TMLE || iterTMLE) {
    IC_dt <- OData$dat.sVar[, list("EIC_i_tplus" = sum(eval(as.name("EIC_i_t")))), by = eval(nodes$IDnode)]
    IC_dt[, ("EIC_i_t0") := res_lastPredQ - mean_est_t]
    IC_dt[, ("EIC_i") := EIC_i_t0 + EIC_i_tplus]
    IC_dt[, c("EIC_i_t0", "EIC_i_tplus") :=  list(NULL, NULL)]
    IC_i_onet <- IC_dt[["EIC_i"]]
    # asymptotic variance (var of the EIC):
    IC_Var <- (1 / (OData$nuniqueIDs)) * sum(IC_dt[["EIC_i"]]^2)
    # variance of the TMLE estimate (scaled by n):
    TMLE_Var <- IC_Var / OData$nuniqueIDs
    # SE of the TMLE
    TMLE_SE <- sqrt(TMLE_Var)
    if (gvars$verbose) print("...empirical mean of the estimated EIC: " %+% mean(IC_dt[["EIC_i"]]))
    if (gvars$verbose) print("...estimated TMLE variance: " %+% TMLE_Var)
    resDF_onet <- cbind(resDF_onet, TMLE_Var = TMLE_Var, TMLE_SE = TMLE_SE)
  }
  return(list(IC_i_onet = IC_i_onet, resDF_onet = resDF_onet))
}
