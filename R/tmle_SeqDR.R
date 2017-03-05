
fitSeqDR <- function(OData,
                        tvals,
                        Qforms,
                        intervened_TRT = NULL,
                        intervened_MONITOR = NULL,
                        rule_name = paste0(c(intervened_TRT, intervened_MONITOR), collapse = ""),
                        models = NULL,
                        estimator = stremrOptions("estimator"),
                        fit_method = stremrOptions("fit_method"),
                        fold_column = stremrOptions("fold_column"),
                        stratifyQ_by_rule = FALSE,
                        useonly_t_TRT = NULL,
                        useonly_t_MONITOR = NULL,
                        CVTMLE = FALSE,
                        trunc_weights = 10^6,
                        parallel = FALSE,
                        return_fW = FALSE,
                        verbose = getOption("stremr.verbose"), ...) {

  stratify_by_last  <- TRUE ## if stratifying we are always stratifying by last treatment only
  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  if (!is.null(models)) assert_that(is.ModelStack(models))
  if (missing(rule_name)) rule_name <- paste0(c(intervened_TRT,intervened_MONITOR), collapse = "")

  # ------------------------------------------------------------------------------------------------
  # **** Evaluate the uncensored and initialize rule followers (everybody is a follower by default)
  # ------------------------------------------------------------------------------------------------
  OData$uncensored <- OData$eval_uncensored()
  OData$follow_rule <- rep.int(TRUE, nrow(OData$dat.sVar)) # (everybody is a follower by default)
  sVar.exprs <- capture.exprs(...)
  models_control <- c(list(models = models), opt_params = list(sVar.exprs))
  models_control[["estimator"]] <- estimator[1L]
  models_control[["fit_method"]] <- fit_method[1L]
  models_control[["fold_column"]] <- fold_column

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
         eval_stabP = FALSE)

  # ------------------------------------------------------------------------------------------
  # Create a back-up of the observed input gstar nodes (created by user in input data):
  # Will add new columns (that were not backed up yet) TO SAME backup data.table
  # ------------------------------------------------------------------------------------------
  OData$backupNodes(c(intervened_TRT,intervened_MONITOR))

  # ------------------------------------------------------------------------------------------------
  # Define the intervention nodes
  # Modify the observed input intervened_NODE in OData$dat.sVar with values from NodeNames for subset_idx
  # ------------------------------------------------------------------------------------------------
  gstar.A <- defineNodeGstarGComp(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, stratifyQ_by_rule, stratify_by_last = stratify_by_last)
  gstar.N <- defineNodeGstarGComp(OData, intervened_MONITOR, nodes$Nnodes, useonly_t_MONITOR, stratifyQ_by_rule, stratify_by_last = stratify_by_last)
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
  est_name <- "SeqDR"
  tmle.run.res <- try(
    if (parallel) {
      mcoptions <- list(preschedule = FALSE)
      '%dopar%' <- foreach::'%dopar%'
      res_byt <- foreach::foreach(t_idx = rev(seq_along(tvals)), .options.multicore = mcoptions) %dopar% {
        t_period <- tvals[t_idx]
        res <- fitSeqDR_onet(OData, t_period, Qforms, stratifyQ_by_rule, CVTMLE = CVTMLE, models = models_control, return_fW = return_fW, verbose = verbose)
        return(res)
      }
      res_byt[] <- res_byt[rev(seq_along(tvals))] # re-assign to order results by increasing t
    } else {
      res_byt <- vector(mode = "list", length = length(tvals))
      for (t_idx in rev(seq_along(tvals))) {
        t_period <- tvals[t_idx]
        cat("Estimating parameter E[Y_d(t)] for t = " %+% t_period, "\n")
        res <- fitSeqDR_onet(OData, t_period, Qforms, stratifyQ_by_rule, CVTMLE = CVTMLE, models = models_control, return_fW = return_fW, verbose = verbose)
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
  ## to extract the EIC estimates by time point (no longer returning those separately, only as part of the estimates data.table)
  # ICs_byt <- resultDT[["IC.St"]]
  # IC.Var.S.d <- t(do.call("cbind", ICs_byt))

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

  res_out <- list(
              # est_name = est_name,
              # periods = tvals,
              # IC.Var.S.d = list(IC.S = IC.Var.S.d),
              # nID = OData$nuniqueIDs,
              # wts_data = { if ((TMLE || iterTMLE) && return_wts) { IPWeights } else { NULL } },
              # rule_name = rule_name,
              # trunc_weights = trunc_weights
              estimates = resultDT
              )

  attr(res_out, "estimator_short") <- est_name
  attr(res_out, "estimator_long") <- est_name
  return(res_out)
}

## ------------------------------------------------------------------------------------------------
## New procedure for sequential double robustness
## ------------------------------------------------------------------------------------------------
fitSeqDR_onet <- function(OData,
                          t_period,
                          Qforms,
                          stratifyQ_by_rule,
                          CVTMLE = FALSE,
                          models,
                          return_fW = FALSE,
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
  # Running fitSeqGcomp_onet might conflict with different Qkplus1
  # ------------------------------------------------------------------------------------------------
  # **** G-COMP: Initiate Qkplus1 - (could be multiple if more than one regimen)
  # That column keeps the tabs on the running Q-fit (SEQ G-COMP)
  # ------------------------------------------------------------------------------------------------
  OData$dat.sVar[, ("EIC_i_t") := 0.0] # set the initial (default values of the t-specific and i-specific EIC estimates)
  OData$dat.sVar[, "Qkplus1" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)
  ## for SDR adding another node:
  # OData$dat.sVar[, "Qstarkprime" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)

  # OData$def.types.sVar() ## was a bottleneck, replaced with below:
  OData$set.sVar.type(name.sVar = "Qkplus1", new.type = "binary")
  # OData$set.sVar.type(name.sVar = "Qstarkprime", new.type = "binary")

  OData$set.sVar.type(name.sVar = "EIC_i_t", new.type = "binary")

  # ------------------------------------------------------------------------------------------------
  # **** Define regression classes for Q.Y and put them in a single list of regressions.
  # **** TO DO: This could also be done only once in the main routine, then just subset the appropriate Q_regs_list
  # ------------------------------------------------------------------------------------------------
  Q_regs_list <- vector(mode = "list", length = length(Qstratas_by_t))
  names(Q_regs_list) <- unlist(Qstratas_by_t)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")

  # SDR_model <- list("objective" = "reg:logistic", "booster" = "gbtree", "nthread" = 1, "max_delta_step" = 6, nrounds = 10)
  SDR_model <- list("objective" = "reg:logistic", "booster" = "gbtree", "nthread" = 1, "max_delta_step" = 6)

  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms_single_t[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    reg <- RegressionClassSDR$new(SDR_model = SDR_model,
                                     Qreg_counter = Qreg_idx[i],
                                     all_Qregs_indx = Qreg_idx,
                                     t_period = Qperiods[i],
                                     TMLE = FALSE, ## set this automatically to FALSE when running SDR:
                                     CVTMLE = CVTMLE,
                                     keep_idx = TRUE, ## Set this automatically to TRUE when running SDR:
                                     stratifyQ_by_rule = stratifyQ_by_rule,
                                     outvar = "Qkplus1",
                                     predvars = regform$predvars,
                                     outvar.class = list("SDRQlearn"), ## Set this automatically to "SDRQlearn" when Running SDR, otherwise "Qlearn"
                                     subset_vars = list("Qkplus1"),
                                     subset_exprs = all_Q_stratify[i],
                                     model_contrl = models,
                                     censoring = FALSE)

    ## For Q-learning this reg class always represents a terminal model class,
    ## since there cannot be any additional model-tree splits by values of subset_vars, subset_exprs, etc.
    ## The following two lines allow for a slightly simplified (shallower) tree representation of GenericModel-type classes.
    ## This also means that stratifying Q fits by some covariate value will not possible with this approach
    ## (i.e., such stratifications would have to be implemented locally by the actual model fitting functions).
    reg_i <- reg$clone()
    reg <- reg_i$ChangeManyToOneRegresssion(1, reg)
    Q_regs_list[[i]] <- reg
  }

  ## TO DO: Automatically call the right constructor below depending on running SDR or regular TMLE
  # Run all Q-learning regressions (one for each subsets defined above, predictions of the last regression form the outcomes for the next:
  Qlearn.fit <- SDRModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)
  # Qlearn.fit$getPsAsW.models()[[1]]
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

  fW_fit <- lastQ.fit$getfit
  resDF_onet[, ("fW_fit") := { if (return_fW) {list(list(fW_fit))} else {list(list(NULL))} }]

  return(resDF_onet)
}