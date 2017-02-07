
## ------------------------------------------------------------------------------------------------
## New procedure for sequential double robustness
## ------------------------------------------------------------------------------------------------
fitSeqDR_onet <- function(OData,
                          t_period,
                          Qforms,
                          stratifyQ_by_rule,
                          models,
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
  names(Qstratas_by_t) <- rep.int("Q.kplus1", length(Qstratas_by_t))
  all_Q_stratify <- Qstratas_by_t

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
  OData$dat.sVar[, ("EIC_i_t") := 0.0] # set the initial (default values of the t-specific and i-specific EIC estimates)
  OData$dat.sVar[, "Q.kplus1" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)
  ## for SDR adding another node:
  OData$dat.sVar[, "Qstarkprime" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)

  # OData$def.types.sVar() ## was a bottleneck, replaced with below:
  OData$set.sVar.type(name.sVar = "Q.kplus1", new.type = "binary")
  OData$set.sVar.type(name.sVar = "EIC_i_t", new.type = "binary")

  # ------------------------------------------------------------------------------------------------
  # **** Define regression classes for Q.Y and put them in a single list of regressions.
  # **** TO DO: This could also be done only once in the main routine, then just subset the appropriate Q_regs_list
  # ------------------------------------------------------------------------------------------------
  Q_regs_list <- vector(mode = "list", length = length(Qstratas_by_t))
  names(Q_regs_list) <- unlist(Qstratas_by_t)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")

  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms_single_t[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    reg <- RegressionClassQlearn$new(Qreg_counter = Qreg_idx[i],
                                     all_Qregs_indx = Qreg_idx,
                                     t_period = Qperiods[i],
                                     TMLE = FALSE, ## set this automatically to FALSE when running SDR:
                                     keep_idx = TRUE, ## Set this automatically to TRUE when running SDR:
                                     stratifyQ_by_rule = stratifyQ_by_rule,
                                     outvar = "Q.kplus1",
                                     predvars = regform$predvars,
                                     outvar.class = list("SDRQlearn"), ## Set this automatically to "SDRQlearn" when Running SDR, otherwise "Qlearn"
                                     subset_vars = list("Q.kplus1"),
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

  browser()

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

  ## 1b. Can instead grab it directly from the data, using the appropriate strata-subsetting expression
  ## mean_est_t and mean_est_t_2 have to be equal!!!!
  # lastQ.fit <- allQmodels[[lastQ_inx]]$getPsAsW.models()[[1]] ## allQmodels[[lastQ_inx]]$get.fits()
  lastQ.fit <- allQmodels[[lastQ_inx]]
  subset_vars <- lastQ.fit$subset_vars
  subset_exprs <- lastQ.fit$subset_exprs
  subset_idx <- OData$evalsubst(subset_vars = subset_vars, subset_exprs = subset_exprs)
  mean_est_t_2 <- mean(OData$dat.sVar[subset_idx, ][["Q.kplus1"]])

  if (gvars$verbose) {
    print("Surv est: " %+% (1 - mean_est_t))
    print("Surv est 2: " %+% (1 - mean_est_t_2))
    print("No. of obs for last prediction of Q: " %+% length(res_lastPredQ))
    print("EY^* estimate at t="%+%t_period %+%": " %+% round(mean_est_t, 5))
  }

  resDF_onet <- data.table(time = t_period,
                           St.GCOMP = NA,
                           St.TMLE = NA,
                           St.iterTMLE = NA,
                           ALLsuccessTMLE = ALLsuccessTMLE,
                           nFailedUpdates = nFailedUpdates,
                           type = ifelse(stratifyQ_by_rule, "stratified", "pooled")
                          )

  est_name <- ifelse(TMLE, "St.TMLE", "St.GCOMP")
  resDF_onet[, (est_name) := (1 - mean_est_t)]

  return(resDF_onet)
}