# ---------------------------------------------------------------------------------------
#' Fit the sequential double robust (SDR) procedure, either with DR tranformation or TMLE-like targeting
#'
#' Interventions on up to 3 nodes are allowed: \code{CENS}, \code{TRT} and \code{MONITOR}.
#' Adjustment will be based on the inverse of the propensity score fits for the observed likelihood (g0.C, g0.A, g0.N),
#' multiplied by the indicator of not being censored and the probability of each intervention in \code{intervened_TRT} and \code{intervened_MONITOR}.
#' Requires column name(s) that specify the counterfactual node values or the counterfactual probabilities of each node being 1 (for stochastic interventions).
#' @param OData Input data object created by \code{importData} function.
#' @param tvals Vector of time-points in the data for which the survival function (and risk) should be estimated
#' @param Qforms Regression formulas, one formula per Q. Only main-terms are allowed.
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
#' @param models Optional parameters specifying the models for fitting the iterative (sequential) G-Computation formula.
#' Must be an object of class \code{ModelStack} specified with \code{gridisl::defModel} function.
#' @param fit_method Model selection approach. Can be either \code{"none"} - no model selection or
#' \code{"cv"} - V fold cross-validation that selects the best model according to lowest cross-validated MSE (must specify the column name that contains the fold IDs).
# \code{"holdout"} - model selection by splitting the data into training and validation samples according to lowest validation sample MSE (must specify the column of \code{TRUE} / \code{FALSE} indicators,
# where \code{TRUE} indicates that this row will be selected as part of the model validation sample).
#' @param fold_column The column name in the input data (ordered factor) that contains the fold IDs to be used as part of the validation sample.
#' Use the provided function \code{\link{define_CVfolds}} to
#' define such folds or define the folds using your own method.
#' @param CVTMLE Set to \code{TRUE} to run the CV-TMLE algorithm instead of the usual TMLE algorithm.
#' Must set either \code{TMLE}=\code{TRUE} or \code{iterTMLE}=\code{TRUE} for this argument to have any effect..
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @param weights Optional \code{data.table} with additional observation- and time-specific weights.  Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}). Not implemented yet.
#' @param parallel Set to \code{TRUE} to run the sequential G-COMP or TMLE in parallel (uses \code{foreach} with \code{dopar} and
#' requires a previously defined parallel back-end cluster)
#' @param return_fW When \code{TRUE}, will return the object fit for the last Q regression as part of the output table.
#' Can be used for obtaining subject-specific predictions of the counterfactual functional E(Y_{d}|W_i).
#' @param use_DR_transform Apply DR transform estimator instead of the iTMLE.
#' @param stabilize Only applies when \code{use_DR_transform=TRUE}. Set this argument to \code{TRUE} to stabilize the weights by the
#' empirical conditional probability of having followed the rule at time-point \code{t}, given the subject has followed the rule all
#' the way up to time-point \code{t}.
#' @param reg_Q (ADVANCED USE ONLY) Directly specify the Q regressions, separately for each time-point.
#' @param SDR_model The xgboost parameter settings for iTMLE non-parametric regression targeting.
#' If missing/NULL the default parameter settings will be used.
#' @param verbose Set to \code{TRUE} to print auxiliary messages during model fitting.
#' @param ... When \code{models} arguments is NOT specified, these additional arguments will be passed on directly to all \code{GridSL}
#' modeling functions that are called from this routine,
#' e.g., \code{family = "binomial"} can be used to specify the model family.
#' Note that all such arguments must be named.
#' @return An output list containing the \code{data.table} with survival estimates over time saved as \code{"estimates"}.
#' @export
fit_iTMLE <- function(OData,
                      tvals,
                      Qforms,
                      intervened_TRT = NULL,
                      intervened_MONITOR = NULL,
                      rule_name = paste0(c(intervened_TRT, intervened_MONITOR), collapse = ""),
                      models = NULL,
                      fit_method = stremrOptions("fit_method"),
                      fold_column = stremrOptions("fold_column"),
                      stratifyQ_by_rule = FALSE,
                      stratify_by_last = TRUE,
                      useonly_t_TRT = NULL,
                      useonly_t_MONITOR = NULL,
                      CVTMLE = FALSE,
                      trunc_weights = 10^6,
                      weights = NULL,
                      parallel = FALSE,
                      return_fW = FALSE,
                      use_DR_transform = FALSE,
                      stabilize = FALSE,
                      reg_Q = NULL,
                      SDR_model = NULL,
                      verbose = getOption("stremr.verbose"), ...) {

  if (is.null(SDR_model)) {
    SDR_model <- list("objective" = "reg:logistic", "booster" = "gbtree", "nthread" = 1, "max_delta_step" = 6, learning_rate = .1, nrounds = 50)
  }

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  assert_that(is.ModelStack(models) || is(models, "Lrnr_base"))
  if (missing(rule_name)) rule_name <- paste0(c(intervened_TRT,intervened_MONITOR), collapse = "")

  # ------------------------------------------------------------------------------------------------
  # **** Evaluate the uncensored and initialize rule followers (everybody is a follower by default)
  # ------------------------------------------------------------------------------------------------
  OData$uncensored <- OData$eval_uncensored()
  OData$follow_rule <- rep.int(TRUE, nrow(OData$dat.sVar)) # (everybody is a follower by default)

  opt_params <- capture.exprs(...)
  if (!("family" %in% names(opt_params))) opt_params[["family"]] <- quasibinomial()

  if (!is.null(models)) {
    assert_that(is.ModelStack(models) || is(models, "Lrnr_base"))
  } else {
    models <- do.call(sl3::Lrnr_glm_fast$new, opt_params)
  }

  models_control <- c(list(models = models), 
                      list(reg_Q = reg_Q), 
                      opt_params = list(opt_params))

  models_control[["fit_method"]] <- fit_method[1L]
  models_control[["fold_column"]] <- fold_column

  if (is.null(OData$fold_column) && is.null(fold_column)) {
    stop("must specify the integer / factor fold_column with validation folds")
  }

  if (!is.null(fold_column)) {
    OData$define_CVfolds(fold_column = fold_column)
  }

  if (missing(tvals)) stop("must specify survival 'tvals' of interest (time period values from column " %+% nodes$tnode %+% ")")

  # ------------------------------------------------------------------------------------------------
  # The weights function is called inside the SDR OBJECT k-LOOP for each [k] row value
  # ------------------------------------------------------------------------------------------------
  OData$IPWeights_info <-
    list(intervened_TRT = intervened_TRT,
         intervened_MONITOR = intervened_MONITOR,
         useonly_t_TRT = useonly_t_TRT,
         useonly_t_MONITOR = useonly_t_MONITOR,
         rule_name = rule_name,
         trunc_weights = trunc_weights,
         holdout = CVTMLE,
         eval_stabP = stabilize)

  # ------------------------------------------------------------------------------------------
  # Create a back-up of the observed input gstar nodes (created by user in input data):
  # Will add new columns (that were not backed up yet) TO SAME backup data.table
  # ------------------------------------------------------------------------------------------
  OData$backupNodes(c(intervened_TRT,intervened_MONITOR))

  # ------------------------------------------------------------------------------------------------
  # Define the intervention nodes
  # Modify the observed input intervened_NODE in OData$dat.sVar with values from NodeNames for subset_idx
  # ------------------------------------------------------------------------------------------------
  gstar.A <- defineNodeGstarGCOMP(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, stratifyQ_by_rule, stratify_by_last = stratify_by_last)
  gstar.N <- defineNodeGstarGCOMP(OData, intervened_MONITOR, nodes$Nnodes, useonly_t_MONITOR, stratifyQ_by_rule, stratify_by_last = stratify_by_last)
  interventionNodes.g0 <- c(nodes$Anodes, nodes$Nnodes)
  interventionNodes.gstar <- c(gstar.A, gstar.N)

  # When the same node names belongs to both g0 and gstar it doesn't need to be intervened upon, so exclude
  common_names <- intersect(interventionNodes.g0, interventionNodes.gstar)
  interventionNodes.g0 <- interventionNodes.g0[!interventionNodes.g0 %in% common_names]
  interventionNodes.gstar <- interventionNodes.gstar[!interventionNodes.gstar %in% common_names]

  OData$interventionNodes.g0 <- interventionNodes.g0
  OData$interventionNodes.gstar <- interventionNodes.gstar

  ## ------------------------------------------------------------------------------------------------
  ## RUN seqDR FOR SEVERAL TIME-POINTS EITHER IN PARALLEL OR SEQUENTIALLY
  ## For est of S(t) over vector of ts, estimate for highest t first going down to smallest t
  ## This is a more efficient when parallelizing, since larger t implies more model runs & longer run time
  ## ------------------------------------------------------------------------------------------------
  est_name <- ifelse(use_DR_transform, "DRtransform", "iTMLE")
  tmle.run.res <- try(
    if (parallel) {
      mcoptions <- list(preschedule = FALSE)
      '%dopar%' <- foreach::'%dopar%'
      res_byt <- foreach::foreach(t_idx = rev(seq_along(tvals)), .options.multicore = mcoptions) %dopar% {
        t_period <- tvals[t_idx]
        res <- fit_iTMLE_onet(OData, t_period, Qforms, stratifyQ_by_rule, CVTMLE = CVTMLE, models = models_control, return_fW = return_fW, use_DR_transform = use_DR_transform, stabilize = stabilize, SDR_model = SDR_model, verbose = verbose)
        return(res)
      }
      res_byt[] <- res_byt[rev(seq_along(tvals))] # re-assign to order results by increasing t
    } else {
      res_byt <- vector(mode = "list", length = length(tvals))
      for (t_idx in rev(seq_along(tvals))) {
        t_period <- tvals[t_idx]
        cat("Estimating parameter E[Y_d(t)] for t = " %+% t_period, "\n")
        res <- fit_iTMLE_onet(OData, t_period, Qforms, stratifyQ_by_rule, CVTMLE = CVTMLE, models = models_control, return_fW = return_fW, use_DR_transform = use_DR_transform, stabilize = stabilize, SDR_model = SDR_model, verbose = verbose)
        res_byt[[t_idx]] <- res
      }
    }
  )

  # ------------------------------------------------------------------------------------------------
  # Restore backed up nodes, even in the event of failure (otherwise the input data is corrupted for good)
  # ------------------------------------------------------------------------------------------------
  OData$restoreNodes(c(intervened_TRT,intervened_MONITOR))

  if (inherits(tmle.run.res, "try-error")) {
    stop(
"...attempt at running TMLE for one or several time-points has failed for unknown reason;
If this error cannot be fixed, consider creating a replicable example and filing a bug report at:
  https://github.com/osofr/stremr/issues
", call. = TRUE)
  }

  resultDT <- rbindlist(res_byt)
  resultDT <- cbind(est_name = est_name, resultDT)
  resultDT[, "rule.name" := eval(as.character(rule_name))]
  attr(resultDT, "estimator_short") <- est_name
  attr(resultDT, "estimator_long") <- est_name
  attr(resultDT, "nID") <- OData$nuniqueIDs
  attr(resultDT, "rule_name") <- rule_name
  attr(resultDT, "stratifyQ_by_rule") <- stratifyQ_by_rule
  attr(resultDT, "stratify_by_last") <- stratify_by_last
  attr(resultDT, "trunc_weights") <- trunc_weights
  attr(resultDT, "time") <- resultDT[["time"]]
  res_out <- list(estimates = resultDT)
  attr(res_out, "estimator_short") <- est_name
  attr(res_out, "estimator_long") <- est_name
  return(res_out)
}

## ------------------------------------------------------------------------------------------------
## New procedure for sequential double robustness
## ------------------------------------------------------------------------------------------------
fit_iTMLE_onet <- function(OData,
                        t_period,
                        Qforms,
                        stratifyQ_by_rule,
                        CVTMLE = FALSE,
                        models,
                        return_fW = FALSE,
                        use_DR_transform = FALSE,
                        stabilize = FALSE,
                        SDR_model,
                        verbose = getOption("stremr.verbose"),
                        ...) {
  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names

  # ------------------------------------------------------------------------------------------------
  # Defining the t periods to loop over FOR A SINGLE RUN OF THE iterative G-COMP/TMLE (one survival point)
  # **** TO DO: The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------
  Qperiods <- rev(OData$min.t:t_period)
  Qreg_idx <- rev(seq_along(Qperiods))
  Qstratas_by_t <- as.list(nodes[['tnode']] %+% " == " %+% (Qperiods))
  names(Qstratas_by_t) <- rep.int("Qkplus1", length(Qstratas_by_t))
  all_Q_stratify <- Qstratas_by_t

  # ------------------------------------------------------------------------------------------------
  # **** Process the input formulas and stratification settings
  # **** TO DO: Add checks that Qforms has correct number of regressions in it
  # ------------------------------------------------------------------------------------------------
  Qforms.default <- rep.int("Qkplus1 ~ Lnodes + Anodes + Cnodes + Nnodes", length(Qperiods))
  if (missing(Qforms)) {
    Qforms_single_t <- Qforms.default
  } else {
    # Qforms_single_t <- Qforms[seq_along(Qperiods)]
    Qforms_single_t <- Qforms[Qreg_idx]
  }

  # ------------------------------------------------------------------------------------------------
  #  ****** Qkplus1 THIS NEEDS TO BE MOVED OUT OF THE MAIN data.table *****
  # Running fit_GCOMP_onet might conflict with different Qkplus1
  # ------------------------------------------------------------------------------------------------
  # **** G-COMP: Initiate Qkplus1 - (could be multiple if more than one regimen)
  # That column keeps the tabs on the running Q-fit (SEQ G-COMP)
  # ------------------------------------------------------------------------------------------------
  OData$dat.sVar[, ("EIC_i_t") := 0.0] # set the initial (default values of the t-specific and i-specific EIC estimates)
  OData$dat.sVar[, ("EIC_i_t_sum") := 0.0] # set the initial rolling EIC sum
  OData$dat.sVar[, "Qkplus1.protected" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)
  OData$dat.sVar[, "Qkplus1" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)
  if ("Qk_hat" %in% names(OData$dat.sVar)) {
    OData$dat.sVar[, "Qk_hat" := NULL]
  }

  # ------------------------------------------------------------------------------------------------
  # **** Define regression classes for Q.Y and put them in a single list of regressions.
  # **** TO DO: This could also be done only once in the main routine, then just subset the appropriate Q_regs_list
  # ------------------------------------------------------------------------------------------------
  Q_regs_list <- vector(mode = "list", length = length(Qstratas_by_t))
  names(Q_regs_list) <- unlist(Qstratas_by_t)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")

  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms_single_t[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    if (!is.null(models[["reg_Q"]])) {
      models[["models"]] <- models[["reg_Q"]][[Qreg_idx[i]]]
    }
    reg <- RegressionClassSDR$new(SDR_model = SDR_model,
                                  stabilize = stabilize,
                                  Qreg_counter = Qreg_idx[i],
                                   all_Qregs_indx = Qreg_idx,
                                   t_period = Qperiods[i],
                                   TMLE = FALSE, ## set this automatically to FALSE when running SDR:
                                   CVTMLE = CVTMLE,
                                   keep_idx = TRUE, ## Set this automatically to TRUE when running SDR:
                                   stratifyQ_by_rule = stratifyQ_by_rule,
                                   outvar = "Qkplus1",
                                   predvars = regform$predvars,
                                   # outvar.class = list("SplitCVSDRQlearn"), ## Set this automatically to "SplitCVSDRQlearn" when Running SDR, otherwise "Qlearn"
                                   # outvar.class = list("SDRtransformQModel"), ## Set this automatically to "SplitCVSDRQlearn" when Running SDR, otherwise "Qlearn"
                                   outvar.class = list(ifelse(use_DR_transform, "SDRtransformQModel", "SplitCVSDRQlearn")),
                                   subset_vars = list("Qkplus1"),
                                   subset_exprs = all_Q_stratify[i],
                                   model_contrl = models,
                                   censoring = FALSE)

    ## For Q-learning this reg class always represents a terminal model class,
    ## since there cannot be any additional model-tree splits by values of subset_vars, subset_exprs, etc.
    ## The following two lines allow for a slightly simplified (shallower) tree representation of ModelGeneric-type classes.
    ## This also means that stratifying Q fits by some covariate value will not possible with this approach
    ## (i.e., such stratifications would have to be implemented locally by the actual model fitting functions).
    reg_i <- reg$clone()
    reg <- reg_i$ChangeManyToOneRegresssion(1, reg)
    Q_regs_list[[i]] <- reg
  }

  ## Automatically call the right constructor below depending on running DR tranform or SDR TMLE
  # Run all Q-learning regressions (one for each subsets defined above, predictions of the last regression form the outcomes for the next:
  if (use_DR_transform) {
    Qlearn.fit <- SDRtransform$new(reg = Q_regs_list, DataStorageClass.g0 = OData)
  } else {
    Qlearn.fit <- SDRModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)
  }

  Qlearn.fit$fit(data = OData)
  OData$Qlearn.fit <- Qlearn.fit

  # get the individual TMLE updates and evaluate if any updates have failed
  allQmodels <- Qlearn.fit$getPsAsW.models()

  allTMLEfits <- lapply(allQmodels, function(Qmod) Qmod$getTMLEfit)
  TMLEfits <- unlist(allTMLEfits)
  successTMLEupdates <- !is.na(TMLEfits) & !is.nan(TMLEfits)
  ALLsuccessTMLE <- all(successTMLEupdates)
  nFailedUpdates <- sum(!successTMLEupdates)

  ## 1a. Grab last reg predictions from Q-regression objects:
  lastQ_inx <- Qreg_idx[1] # The index for the last Q-fit (first time-point)

  ## Get the previously saved mean prediction for Q from the very last regression (first time-point, all n obs):
  res_lastPredQ <- allQmodels[[lastQ_inx]]$predictAeqa()
  mean_est_t <- mean(res_lastPredQ)

  ## 1b. Can instead grab the last prediction (Qk_hat) directly from the data, using the appropriate strata-subsetting expression
  ## mean_est_t and mean_est_t_2 have to be equal!!!!
  # lastQ.fit <- allQmodels[[lastQ_inx]]$getPsAsW.models()[[1]] ## allQmodels[[lastQ_inx]]$get.fits()
  lastQ.fit <- allQmodels[[lastQ_inx]]
  subset_vars <- lastQ.fit$subset_vars
  subset_exprs <- lastQ.fit$subset_exprs
  subset_idx <- OData$evalsubst(subset_vars = subset_vars, subset_exprs = subset_exprs)
  mean_est_t_2 <- mean(OData$dat.sVar[subset_idx, ][["Qk_hat"]])

  if (gvars$verbose) {
    print("Surv est: " %+% (1 - mean_est_t))
    print("Surv est 2: " %+% (1 - mean_est_t_2))
    print("No. of obs for last prediction of Q: " %+% length(res_lastPredQ))
    print("EY^* estimate at t="%+%t_period %+%": " %+% round(mean_est_t, 5))
  }

  resDF_onet <- data.table(time = t_period,
                           St.SDR = NA,
                           type = ifelse(stratifyQ_by_rule, "stratified", "pooled")
                          )

  est_name <- "St.SDR"
  resDF_onet[, (est_name) := (1 - mean_est_t)]

  # ------------------------------------------------------------------------------------------------
  # SDR INFERENCE
  # ------------------------------------------------------------------------------------------------
  IC_i_onet <- vector(mode = "numeric", length = OData$nuniqueIDs)
  IC_i_onet[] <- NA
  IC_dt <- OData$dat.sVar[, list("EIC_i_tplus" = sum(eval(as.name("EIC_i_t")))), by = eval(nodes$IDnode)]
  IC_dt[, ("EIC_i_t0") := res_lastPredQ - mean_est_t]

  ## FOR DR-transform the i-specific EIC estimate (over all t) is just (res_lastPredQ - psi_n)
  ## This is because the very last Qhat (res_lastPredQ) is based on the transformed Gamma(k) for k=0 (DR transform for Q_0),
  ## which over-wrote the initial Q_hat at first time-point.
  ## The Gamma(0) was evaluated as sum_{k}[wt(k)(Q_{k+1}-Q_{k})] + Q_hat.
  ## In the last sum Q_hat was the initial regression fit E[Gamma(1)|h(0),A(0)=1] for the very first time-point.
  ## Thus, but subtracting the mean parameter estimate (psi_n) from res_lastPredQ we get the i-specific EIC estimate.
  if (use_DR_transform) {
    IC_dt[, ("EIC_i") := EIC_i_t0]
  } else {
    IC_dt[, ("EIC_i") := EIC_i_t0 + EIC_i_tplus]
  }

  IC_dt[, c("EIC_i_t0", "EIC_i_tplus") :=  list(NULL, NULL)]
  IC_i_onet <- IC_dt[["EIC_i"]]
  ## estimate of the asymptotic variance (var of the EIC):
  IC_Var <- (1 / (OData$nuniqueIDs)) * sum(IC_dt[["EIC_i"]]^2)
  ## estimate of the variance of SDR estimate (scaled by n):
  SDR_Var <- IC_Var / OData$nuniqueIDs
  ## SE of the SDR
  SDR_SE <- sqrt(SDR_Var)
  if (gvars$verbose) {
    print("...empirical mean of the estimated EIC: " %+% mean(IC_dt[["EIC_i"]]))
    print("...estimated SDR variance: " %+% SDR_Var)
  }
  resDF_onet[, ("SE.SDR") := SDR_SE]
  ## save the i-specific estimates of the EIC as a separate column:
  resDF_onet[, ("IC.St") := list(list(IC_i_onet))]
  fW_fit <- lastQ.fit$getfit
  resDF_onet[, ("fW_fit") := { if (return_fW) {list(list(fW_fit))} else {list(list(NULL))} }]
  return(resDF_onet)
}