fit_GCOMP_onet <- function(OData,
                             t_period,
                             Qforms,
                             Qstratify,
                             stratifyQ_by_rule,
                             TMLE,
                             iterTMLE,
                             CVTMLE = FALSE,
                             byfold_Q = FALSE,
                             models,
                             max_iter = 10,
                             adapt_stop = TRUE,
                             adapt_stop_factor = 10,
                             tol_eps = 0.001,
                             return_fW = FALSE,
                             maxpY = 1.0,
                             TMLE_updater = "TMLE.updater.speedglm",
                             verbose = getOption("stremr.verbose")) {
  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names

  # ------------------------------------------------------------------------------------------------
  # Defining the t periods to loop over FOR A SINGLE RUN OF THE iterative G-COMP/TMLE (one survival point)
  # **** TO DO ****: The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------
  Qperiods <- rev(OData$min.t:t_period)
  Qreg_idx <- rev(seq_along(Qperiods))
  Qstratas_by_t <- as.list(nodes[['tnode']] %+% " == " %+% (Qperiods))
  names(Qstratas_by_t) <- rep.int("Qkplus1", length(Qstratas_by_t))

  # Adding user-specified stratas to each Q(t) regression:
  all_Q_stratify <- Qstratas_by_t
  if (!is.null(Qstratify)) {
    stop("...Qstratify is not implemented...")
    assert_that(is.vector(Qstratify))
    assert_that(is.character(Qstratify))
    for (idx in seq_along(all_Q_stratify)) {
      all_Q_stratify[[idx]] <- paste0(all_Q_stratify[[idx]], " & ", Qstratify)
    }
  }

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
  OData$dat.sVar[, "Qkplus1" := as.numeric(get(OData$nodes$Ynode))] # set the initial values of Q (the observed outcome node)

  if ("Qk_hat" %in% names(OData$dat.sVar)) {
    OData$dat.sVar[, "Qk_hat" := NULL]
  }
  if ("res_lastPredQ" %in% names(OData$dat.sVar)) {
    OData$dat.sVar[, "res_lastPredQ" := NULL]
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
    last.reg <- Qreg_idx[i] == 1L
    reg <- RegressionClassQlearn$new(Qreg_counter = Qreg_idx[i],
                                     all_Qregs_indx = Qreg_idx,
                                     t_period = Qperiods[i],
                                     TMLE = TMLE,
                                     CVTMLE = CVTMLE,
                                     byfold_Q = byfold_Q,
                                     keep_idx = ifelse(iterTMLE, TRUE, FALSE),
                                     # keep_model_fit = ifelse(iterTMLE, TRUE, last.reg),
                                     keep_model_fit = ifelse(iterTMLE, TRUE, FALSE),
                                     stratifyQ_by_rule = stratifyQ_by_rule,
                                     outvar = "Qkplus1",
                                     predvars = regform$predvars,
                                     outvar.class = list("Qlearn"),
                                     subset_vars = list("Qkplus1"),
                                     subset_exprs = all_Q_stratify[i],
                                     model_contrl = models,
                                     censoring = FALSE,
                                     maxpY = maxpY,
                                     TMLE_updater = TMLE_updater)

    ## For Q-learning this reg class always represents a terminal model class,
    ## since there cannot be any additional model-tree splits by values of subset_vars, subset_exprs, etc.
    ## The following two lines allow for a slightly simplified (shallower) tree representation of ModelGeneric-type classes.
    ## This also means that stratifying Q fits by some covariate value will not possible with this approach
    ## (i.e., such stratifications would have to be implemented locally by the actual model fitting functions).
    reg_i <- reg$clone()
    reg <- reg_i$ChangeManyToOneRegresssion(1, reg)
    Q_regs_list[[i]] <- reg
  }

  # Run all Q-learning regressions (one for each subsets defined above, predictions of the last regression form the outcomes for the next:
  Qlearn.fit <- ModelGeneric$new(reg = Q_regs_list, DataStorageClass.g0 = OData)
  Qlearn.fit$fit(data = OData, Qlearn.fit = Qlearn.fit)
  # Qlearn.fit$getPsAsW.models()[[1]]
  OData$Qlearn.fit <- Qlearn.fit

  # When the model is fit with user-defined stratas, need special functions to extract the final fit (current aproach will not work):
  # allQmodels[[1]]$getPsAsW.models()[[1]]$getPsAsW.models()

  # get the individual TMLE updates and evaluate if any updates have failed
  allQmodels <- Qlearn.fit$getPsAsW.models()

  lastQ_inx <- Qreg_idx[1] # The index for the last Q-fit (first time-point)
  lastQ.fit <- allQmodels[[lastQ_inx]]
  subset_vars <- lastQ.fit$subset_vars
  subset_exprs <- lastQ.fit$subset_exprs
  subset_idx <- OData$evalsubst(subset_exprs = subset_exprs) ## subset_vars = subset_vars,

  ## Get the previously saved mean prediction for Q from the very last regression (first time-point, all n obs):
  res_lastPredQ <- lastQ.fit$predictAeqa()
  ## remove vector of predictions from the modeling class (to free up memory):
  lastQ.fit$wipe.probs

  ## Can instead grab the last prediction of Q (Qk_hat) directly from the data, using the appropriate strata-subsetting expression:
  mean_est_t_2 <- OData$dat.sVar[subset_idx, list("cum.inc" = mean(Qk_hat)), by = "rule.name"]

  ## Evaluate the last predictions by rule (which rule each prediction belongs to):
  OData$dat.sVar[subset_idx, "res_lastPredQ" := res_lastPredQ]
  mean_est_t <- OData$dat.sVar[subset_idx, list("cum.inc" = mean(res_lastPredQ)), by = "rule.name"]

  if (gvars$verbose) {
    print("mean est: ");  print(mean_est_t)
    print("mean est 2: "); print(mean_est_t_2)
    print("No. of obs for last prediction of Q: " %+% length(res_lastPredQ))
  }

  resDF_onet <- data.table(time = t_period,
                           St.GCOMP = NA,
                           St.TMLE = NA,
                           type = ifelse(stratifyQ_by_rule, "stratified", "pooled"),
                           rule.name = unique(OData$dat.sVar[["rule.name"]])
                          )

  est_name <- ifelse(TMLE, "St.TMLE", "St.GCOMP")
  resDF_onet <- resDF_onet[mean_est_t_2, on = "rule.name"]

  resDF_onet[, (est_name) := (1 - cum.inc)]

  # ------------------------------------------------------------------------------------------------
  # RUN ITERATIVE TMLE (updating all Q's at once):
  # ------------------------------------------------------------------------------------------------
  if (iterTMLE){
    res <- iterTMLE_onet(OData, Qlearn.fit, Qreg_idx, max_iter = max_iter, adapt_stop = adapt_stop, adapt_stop_factor = adapt_stop_factor, tol_eps = tol_eps)
    ## 1a. Grab the mean prediction from the very last regression (over all n observations);
    lastQ_inx <- Qreg_idx[1]
    res_lastPredQ <- allQmodels[[lastQ_inx]]$predictAeqa()

    mean_est_t <- mean(res_lastPredQ)
    if (gvars$verbose) print("Iterative TMLE surv estimate: " %+% (1 - mean_est_t))

    ## 1b. Grab it directly from the data, using the appropriate strata-subsetting expression
    resDF_onet[, ("St.TMLE") := (1 - mean_est_t)]
  }

  # ------------------------------------------------------------------------------------------------
  # TMLE INFERENCE
  # ------------------------------------------------------------------------------------------------
  IC_i_onet <- vector(mode = "numeric", length = OData$nuniqueIDs)
  IC_i_onet[] <- NA
  ## save the i-specific estimates of the EIC as a separate column:
  resDF_onet[, ("IC.St") := list(list(IC_i_onet))]

  if (TMLE || iterTMLE) {
    resDF_onet[, ("IC.St") := NULL]
    IC_dt <- OData$dat.sVar[, list("EIC_i_tplus" = sum(eval(as.name("EIC_i_t")))), by = c(nodes$IDnode, "rule.name")]

    IC_dt_t0 <- OData$dat.sVar[subset_idx,  c(nodes$IDnode, "res_lastPredQ", "rule.name"), with = FALSE][,
      "mean_Q" := mean(res_lastPredQ), by = c("rule.name")][,
      "EIC_i_t0" := res_lastPredQ-mean_Q][,
      c("mean_Q", "res_lastPredQ") := list(NULL, NULL)
      ]

    ## remove newly added column from observed dataset
    OData$dat.sVar[, "res_lastPredQ" := NULL]

    IC_dt <- IC_dt[IC_dt_t0, on = c(nodes$IDnode, "rule.name")]
    IC_dt[, ("EIC_i") := EIC_i_t0 + EIC_i_tplus]
    IC_dt[, c("EIC_i_t0", "EIC_i_tplus") :=  list(NULL, NULL)]

    ## asymptotic variance (var of the EIC), by rule:
    IC_Var <- IC_dt[,
      list("IC_Var" = (1 / (.N)) * sum(EIC_i^2),
           "TMLE_Var" = (1 / (.N)) * sum(EIC_i^2) / .N
            ), by = "rule.name"][,
      "SE.TMLE" := sqrt(TMLE_Var)]
    IC_Var[, "IC_Var" := NULL][, "TMLE_Var" := NULL]

    if (gvars$verbose) {
      print("...empirical mean of the estimated EIC: ")
      print(IC_dt[, list("mean_EIC" = mean(EIC_i)), by = "rule.name"])
      print("...estimated TMLE variance: "); print(IC_Var)
    }

    resDF_onet <- resDF_onet[IC_Var, on = "rule.name"]
    data.table::setkeyv(IC_dt, c("rule.name", nodes$IDnode))
    IC_i_onet_byrule <- IC_dt[, (nodes$IDnode) := NULL] %>%
                        tibble::as_tibble() %>%
                        tidyr::nest(EIC_i, .key = "IC.St") %>%
                        dplyr::mutate(IC.St = map(IC.St, ~ .x[[1]]))
    data.table::setDT(IC_i_onet_byrule)
    # IC_i_onet_byrule[, "IC.St" := list(list(as.numeric(IC.St[[1]][[1]]))), by = "rule.name"]
    resDF_onet <- resDF_onet[IC_i_onet_byrule, on = "rule.name"]
  }

  fW_fit <- lastQ.fit$getfit
  resDF_onet[, ("fW_fit") := { if (return_fW) {list(list(fW_fit))} else {list(list(NULL))} }]

  ## re-order the columns so that "rule.name" col is always last
  data.table::setcolorder(resDF_onet, c(setdiff(names(resDF_onet), "rule.name"), "rule.name"))
  return(resDF_onet)
}