# ------------------------------------------------------------------------------------------
# H2O ISSUES WITH PARALLEL TRAINING:
# ------------------------------------------------------------------------------------------
# 1) EVERYTIME a subset [,] is done on frame (H2O.dat.sVar[rows_subset, vars]) a new frame is created with an automatically ID
# 2) EVERYTIME a prediction is made a new h2o frame is created: predictions_"modelIDname"_on_"H2OFRAMENAME"
# 3) EVERYTIME a prediction result is pulled ([,"p1"]), a new temp FRAME  is created
# 4) WHEN DOING "subsetH2Oframe[, outvar]" in BinomialH2O$setdataH2O A NEW temp FRAME (AUTOMATIC ID) is CREATED
# 5) Having Q.kplus1 as a column in the main DT or FRAME is also an issue, since parallel training might overwrite it

# ------------------------------------------------------------------------------------------
# POTENTIAL SOLUTION:
# ------------------------------------------------------------------------------------------
# *) PASS Q.kplus1 AS a vector argument?
# *) using existing approach to LOADING A NEW FRAME FOR EVERY single call to QlearnModel$fit,
#    but pass the args so that destination_frame is always different for each separate call of fitSeqGcomp_singlet()
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
# *** Pooling & performing only one TMLE update across all Q.k ***
#     In general, we can pool all Q.kplus (across all t's) and do a single update on all of them
#     Would make this TMLE iterative, but would not require refitting on the initial Q's
#     See how fast is the G-comp formula prediction without fitting (the subseting/prediction loop alone) -> maybe fast enough to allow iterations
# ------------------------------------------------------------------------------------------
# *** Stochastic interventions ***
#     Deal with stochastic interventions in QlearnModel class (possibly daughter classes) with direct intergration or MC sim
#     Need to do either MC integration (sample from g.star then predict Q)
#     or direct weighted sum of predicted Q's with weights given by g.star (on A and N).
# *) gstar_TRT/gstar_MONITOR are vectors of coutnerfactual probs -> allow each to be multivariate (>1 cols)
# *) For column 0<gstar_TRT<1 it defines the counterfactual probability that P(TRT[t]^*=1)=gstar_TRT[t].
# *) gstar can be a vector, i.e., abar=(0,0,0,0) or a matrix of counterfactual treatments, like in ltmle.
# ------------------------------------------------------------------------------------------
#  *** rule_followers:
#     Once you go off the treatment first time, this is it, the person is censored for the rest of the follow-up.
#     Rule followers are now evaluated automatically by comparing (gstar_TRT and TRT and CENS) and (gstar_MONITOR and MONITOR and CENS).
#     Rule followers are: (gstar_TRT = 1) & (TRT == 1) or (gstar_TRT = 0) & (TRT == 0) or (gstar_TRT > 0 & gstar_TRT < 0) & (Not Censored).
#     Exactly the same logic is also applied to (gstar_MONITOR & MONITOR) when these are specified.
#     If either TRT or MONITOR is multivariate (>1 col), this logic needs to be applied FOR EACH COLUMN of TRT/MONITOR.
# ------------------------------------------------------------------------------------------
# *** Accessing correct QlearnModel ***
#     Need to come up with a good way of accessing correct terminal QlearnModel to get final Q-predictions in private$probAeqa.
#     These predictions are for the final n observations that were used for evaluating E[Y^a].
#     However, if we do stratify=TRUE and use a stack version of g-comp with different strata all in one stacked database it will be difficult
#     to pull the right QlearnModel.
#     Alaternative is to ignore that completely and go directly for the observed data (by looking up the right column/row combos).
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
#' @export
fitTMLE <- function(...) {
  fitSeqGcomp(TMLE = TRUE, ...)
}

#' @export
fitSeqGcomp <- function(OData,
                        t_periods,
                        Qforms,
                        gstar_TRT = NULL,
                        gstar_MONITOR = NULL,
                        stratifyQ_by_rule = FALSE,
                        TMLE = FALSE,
                        rule_name = paste0(c(gstar_TRT, gstar_MONITOR), collapse = ""),
                        IPWeights,
                        params_Q = list(),
                        parallel = FALSE,
                        verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  assert_that(is.list(params_Q))

  if (missing(rule_name)) rule_name <- paste0(c(gstar_TRT,gstar_MONITOR), collapse = "")
  # ------------------------------------------------------------------------------------------------
  # **** Evaluate the uncensored and initialize rule followers (everybody is a follower by default)
  # **** NOTE: THIS NEEDS TO BE TAKEN OUT OF HERE AND PUT AS A SEPARATE FUNCTION TO BE CALLED FROM getIPWeights and/or fitSeqGcomp
  # ------------------------------------------------------------------------------------------------
  OData$uncensored_idx <- OData$eval_uncensored()
  OData$rule_followers_idx <- rep.int(TRUE, nrow(OData$dat.sVar)) # (everybody is a follower by default)

  # ------------------------------------------------------------------------------------------------
  # **** Define the intervention nodes
  # ------------------------------------------------------------------------------------------------
  if (!is.null(gstar_TRT)) {
    gstar.A <- gstar_TRT
    for (gstar_TRT_col in gstar_TRT) CheckVarNameExists(OData$dat.sVar, gstar_TRT_col)
    # UPDATE RULE FOLLOWERS FOR TRT IF DOING stratified G-COMP:
    if (stratifyQ_by_rule) {
      rule_followers_idx <- OData$eval_rule_followers(NodeName = nodes$Anodes, gstar.NodeName = gstar.A)
      OData$rule_followers_idx <- rule_followers_idx & OData$rule_followers_idx & OData$uncensored_idx
    }
  } else {
    gstar.A <- nodes$Anodes # use the actual observed exposure (no intervention on TRT)
  }

  if (!is.null(gstar_MONITOR)) {
    gstar.N <- gstar_MONITOR
    for (gstar_MONITOR_col in gstar_MONITOR) CheckVarNameExists(OData$dat.sVar, gstar_MONITOR_col)
    # UPDATE RULE FOLLOWERS FOR MONITOR IF DOING stratified G-COMP:
    if (stratifyQ_by_rule) {
      rule_followers_idx <- OData$eval_rule_followers(NodeName = nodes$Nnodes, gstar.NodeName = gstar.N)
      OData$rule_followers_idx <- rule_followers_idx & OData$rule_followers_idx & OData$uncensored_idx
    }
  } else {
    gstar.N <- nodes$Nnodes # use the actual observed monitoring probability (no intervention on MONITOR)
  }

  interventionNodes.g0 <- c(nodes$Anodes, nodes$Nnodes)
  interventionNodes.gstar <- c(gstar.A, gstar.N)
  OData$interventionNodes.g0 <- interventionNodes.g0
  OData$interventionNodes.gstar <- interventionNodes.gstar




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
  params_Q$solver <- "IRLSM"
  # parameter for running GBM with continuous outcome (default is classification with "bernoulli" and 0/1 outcome):
  params_Q$distribution <- "gaussian"

  # ------------------------------------------------------------------------------------------------
  # **** Add weights if TMLE=TRUE and if weights were defined
  # NOTE: This needs to be done only once if evaluating survival over several t_periods
  # ------------------------------------------------------------------------------------------------
  if (TMLE) {
    if (missing(IPWeights)) {
      message("...evaluating IPWeights for TMLE...")
      IPWeights <- getIPWeights(OData, gstar_TRT, gstar_MONITOR, rule_name, stabilize = FALSE)
    } else {
      getIPWeights_fun_call <- attributes(IPWeights)[['getIPWeights_fun_call']]
      message("applying user-specified IPWeights, make sure these weights were obtained by making a call: \n'getIPWeights((OData, gstar_TRT, gstar_MONITOR, stabilize = FALSE)'")
      message("the currently supplied weights were obtained with a call: \n" %+% deparse(getIPWeights_fun_call)[[1]])
      assert_that(all.equal(attributes(IPWeights)[['gstar_TRT']], gstar_TRT))
      assert_that(all.equal(attributes(IPWeights)[['gstar_MONITOR']], gstar_MONITOR))
      assert_that(attributes(IPWeights)[['stabilize']] == FALSE)
    }
    assert_that(is.data.table(IPWeights))
    assert_that("cumm.IPAW" %in% names(IPWeights))
    OData$IPwts_by_regimen <- IPWeights
  }

  if (missing(t_periods)) stop("must specify survival 't_periods' of interest (time period values from column " %+% nodes$tnode %+% ")")

  # ------------------------------------------------------------------------------------------------
  # RUN GCOMP OR TMLE FOR SEVERAL TIME-POINTS EITHER IN PARALLEL OR SEQUENTIALLY
  # ------------------------------------------------------------------------------------------------
  est_name <- ifelse(TMLE, "TMLE", "GCOMP")
  if (parallel) {
    mcoptions <- list(preschedule = FALSE)
    res_byt <- foreach::foreach(t_idx = seq_along(t_periods), .options.multicore = mcoptions) %dopar% {
                                  t_period <- t_periods[t_idx]
                                  riskP1_byt <- fitSeqGcomp_singlet(OData, t_period, Qforms, stratifyQ_by_rule, TMLE, params_Q, verbose)
                                  surv_byt <- 1-riskP1_byt
                                  data.frame(t = t_period, risk = riskP1_byt, surv = surv_byt)
                }
  } else {
    res_byt <- vector(mode = "list", length = length(t_periods))
    for (t_idx in seq_along(t_periods)) {
      t_period <- t_periods[t_idx]
      riskP1_byt <- fitSeqGcomp_singlet(OData, t_period, Qforms, stratifyQ_by_rule, TMLE, params_Q, verbose)
      surv_byt <- 1-riskP1_byt
      res_byt[[t_idx]] <- data.frame(t = t_period, risk = riskP1_byt, surv = surv_byt)
    }
  }
  resultDT <- data.table(est_name = est_name, rbindlist(res_byt))
  return(resultDT)
}

fitSeqGcomp_singlet <- function(OData,
                                t_period,
                                Qforms,
                                stratifyQ_by_rule,
                                TMLE,
                                params_Q,
                                verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names

  # ------------------------------------------------------------------------------------------------
  # Defining the t periods to loop over FOR A SINGLE RUN OF THE iterative G-COMP/TMLE (one survival point)
  # **** TO DO: The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------
  Qperiods <- rev(OData$min.t:t_period)
  Qreg_idx <- rev(seq_along(Qperiods))
  stratify_Q <- as.list(nodes[['tnode']] %+% " == " %+% (Qperiods))
  names(stratify_Q) <- rep.int("Q.kplus1", length(stratify_Q))

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
  # Running fitSeqGcomp_singlet might conflict with different Q.kplus1
  # ------------------------------------------------------------------------------------------------
  # **** G-COMP: Initiate Q.kplus1 - (could be multiple if more than one regimen)
  # That column keeps the tabs on the running Q-fit (SEQ G-COMP)
  # ------------------------------------------------------------------------------------------------
  OData$dat.sVar[, "Q.kplus1" := as.numeric(get(OData$nodes$Ynode))]
  OData$def.types.sVar()

  # ------------------------------------------------------------------------------------------------
  # **** Define regression classes for Q.Y and put them in a single list of regressions.
  # **** TO DO: This could also be done only once in the main routine, then just subset the appropriate Q_regs_list
  # ------------------------------------------------------------------------------------------------
  Q_regs_list <- vector(mode = "list", length = length(stratify_Q))
  names(Q_regs_list) <- unlist(stratify_Q)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")
  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms_single_t[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    reg <- RegressionClassQlearn$new(Qreg_counter = Qreg_idx[i], t_period = Qperiods[i], stratifyQ_by_rule = stratifyQ_by_rule,
                                     TMLE = TMLE,
                                     outvar = "Q.kplus1", predvars = regform$predvars, outvar.class = list("Qlearn"),
                                     subset_vars = list("Q.kplus1"), subset_exprs = stratify_Q[i], model_contrl = params_Q,
                                     censoring = FALSE)
    Q_regs_list[[i]] <- reg
  }
  Qlearn.fit <- GenericModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)
  # Run all Q-learning regressions (one for each subsets defined above, predictions of the last regression form the outcomes for the next:
  Qlearn.fit$fit(data = OData)
  OData$Qlearn.fit <- Qlearn.fit

  # 1a. Grab the mean prediction from the very last regression (over all n observations);
  # this is the G-COMP/TMLE estimate of survival for a single t period (t specified in the first Q-reg)
  lastQ_inx <- Qreg_idx[1] # the index for the last Q-fit
  res_lastPredQ_Prob1 <- Qlearn.fit$predictRegK(lastQ_inx, OData$nuniqueIDs)
  mean_est_t <- mean(res_lastPredQ_Prob1)
  print("No. of obs for last prediction of Q: " %+% length(res_lastPredQ_Prob1))
  print("EY^* estimate at t="%+%t_period %+%": " %+% round(mean_est_t, 5))
  # # 1b. Grab the right model (QlearnModel) and pull it directly:
  #   lastQ.fit <- Qlearn.fit$getPsAsW.models()[[lastQ_inx]]$getPsAsW.models()[[1]]
  #   lastQ.fit
  #   # for all observations in long format:
  #   length(lastQ.fit$getprobA1)
  #   head(lastQ.fit$getprobA1)
  #   # for all unique IDs in the data:
  #   length(lastQ.fit$predictAeqa())
  #   head(lastQ.fit$predictAeqa())
  #   mean(lastQ.fit$predictAeqa())
  # # 1c. Grab it directly from the data, using the appropriate strata-subsetting expression
  #   subset_vars <- lastQ.fit$subset_vars
  #   subset_exprs <- lastQ.fit$subset_exprs
  #   subset_idx <- OData$evalsubst(subset_vars = subset_vars, subset_exprs = subset_exprs)
  #   mean(OData$dat.sVar[subset_idx, ][["Q.kplus1"]])

  return(mean_est_t)
}


