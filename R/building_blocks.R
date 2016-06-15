
# ------------------------------------------------------------------
# - BLOCK 1: Process inputs and define OData R6 object
# ------------------------------------------------------------------
#' @export
get_Odata <- function(data, ID = "Subj_ID", t.name = "time_period", covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                      noCENS.cat = 0L, verbose = getOption("stremr.verbose")) {
# SHIFTUPoutcome = TRUE,

  gvars$verbose <- verbose
  gvars$noCENS.cat <- noCENS.cat
  if (verbose) {
    current.options <- capture.output(str(gvars$opts))
    print("Using the following stremr options/settings: ")
    cat('\n')
    cat(paste0(current.options, collapse = '\n'), '\n')
  }
  if (missing(covars)) { # define time-varing covars (L) as everything else in data besides these vars
    covars <- setdiff(colnames(data), c(ID, t.name, CENS, TRT, MONITOR, OUTCOME))
  }
  # The ordering of variables in this list is the assumed temporal order!
  nodes <- list(Lnodes = covars, Cnodes = CENS, Anodes = TRT, Nnodes = MONITOR, Ynode = OUTCOME, IDnode = ID, tnode = t.name)
  OData <- DataStorageClass$new(Odata = data, nodes = nodes, noCENS.cat = noCENS.cat)

  # --------------------------------------------------------------------------------------------------------
  # Create dummies for factors
  # --------------------------------------------------------------------------------------------------------
  factor.Ls <- unlist(lapply(OData$dat.sVar, is.factor))
  factor.Ls <- names(factor.Ls)[factor.Ls]
  new.factor.names <- vector(mode="list", length=length(factor.Ls))
  names(new.factor.names) <- factor.Ls
  if (length(factor.Ls)>0 && verbose)
    message("...converting the following factor(s) to binary dummies (and droping the first factor levels): " %+% paste0(factor.Ls, collapse=","))
  for (factor.varnm in factor.Ls) {
    factor.levs <- levels(OData$dat.sVar[,factor.varnm, with=FALSE][[1]])
    factor.levs <- factor.levs[-1] # remove the first level (reference class)
    # use levels to define cat indicators:
    OData$dat.sVar[,(factor.varnm %+% "_" %+% factor.levs) := lapply(factor.levs, function(x) levels(get(factor.varnm))[get(factor.varnm)] %in% x)]
    # to remove the origional factor var: # OData$dat.sVar[,(factor.varnm):=NULL]
    new.factor.names[[factor.varnm]] <- factor.varnm %+% "_" %+% factor.levs
    # alternative wth dcast: # out <- dcast(OData$dat.sVar, "StudyID + intnum + race ~ race", fun = length, value.var = "race")
  }
  OData$new.factor.names <- new.factor.names

  # --------------------------------------------------------------------------------------------------------
  # Convert all logical vars to binary integers
  # --------------------------------------------------------------------------------------------------------
  if (verbose) message("...converting logical columns to binary integers (0 = FALSE)...")
  logical.Ls <- unlist(lapply(OData$dat.sVar, is.logical))
  logical.Ls <- names(logical.Ls)[logical.Ls]
  for (logical.varnm in logical.Ls) {
    OData$dat.sVar[,(logical.varnm) := as.integer(get(logical.varnm))]
  }

  # # ---------------------------------------------------------------------------
  # # DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
  # # Might expand this in the future to allow defining arbitrary summaries
  # # ---------------------------------------------------------------------------
  # lagnodes <- c(nodes$Cnodes, nodes$Anodes, nodes$Nnodes)
  # newVarnames <- lagnodes %+% ".tminus1"
  # OData$dat.sVar[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=get(nodes$ID), .SDcols=(lagnodes)]
  # # -------------------------------------------------------------------------------------------
  # # Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
  # # -------------------------------------------------------------------------------------------
  # Ynode <- nodes$Ynode
  # shifted.OUTCOME <- Ynode%+%".tplus1"
  # if (SHIFTUPoutcome) {
  #   OData$dat.sVar[, (shifted.OUTCOME) := shift(get(Ynode), n = 1L, type = "lead"), by = eval(nodes$IDnode)]
  # } else {
  #   OData$dat.sVar[, (shifted.OUTCOME) := get(Ynode)]
  # }
  # # OData$dat.sVar <- OData$dat.sVar[!is.na(get(shifted.OUTCOME)), ] # drop and over-write previous data.table, removing last rows.

  for (Cnode in nodes$Cnodes) CheckVarNameExists(OData$dat.sVar, Cnode)
  for (Anode in nodes$Anodes) CheckVarNameExists(OData$dat.sVar, Anode)
  for (Nnode in nodes$Nnodes) CheckVarNameExists(OData$dat.sVar, Nnode)
  for (Ynode in nodes$Ynode)  CheckVarNameExists(OData$dat.sVar, Ynode)
  for (Lnode in nodes$Lnodes) CheckVarNameExists(OData$dat.sVar, Lnode)

  return(OData)
}

# ------------------------------------------------------------------
# - BLOCK 2: define regression models, define a single RegressionClass & fit the propensity score for observed data, summary.g0 g0 (C,A,N)
# ------------------------------------------------------------------
#' @export
get_fits <- function(OData, gform.CENS, gform.TRT, gform.MONITOR,
                    stratify.CENS = NULL, stratify.TRT = NULL, stratify.MONITOR = NULL,
                    verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  # ------------------------------------------------------------------------------------------------
  # Process the input formulas and stratification settings;
  # Define regression classes for g.C, g.A, g.N and put them in a single list of regressions.
  # ------------------------------------------------------------------------------------------------
  gform.CENS.default <- "Cnodes ~ Lnodes"
  gform.TRT.default <- "Anodes ~ Lnodes"
  gform.MONITOR.default <- "Nnodes ~ Anodes + Lnodes"
  g_CAN_regs_list <- vector(mode = "list", length = 3)
  names(g_CAN_regs_list) <- c("gC", "gA", "gN")
  class(g_CAN_regs_list) <- c(class(g_CAN_regs_list), "ListOfRegressionForms")
  gC.sVars <- process_regforms(regforms = gform.CENS, default.reg = gform.CENS.default, stratify.EXPRS = stratify.CENS,
                              OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = TRUE)
  g_CAN_regs_list[["gC"]] <- gC.sVars$regs

  gA.sVars <- process_regforms(regforms = gform.TRT, default.reg = gform.TRT.default, stratify.EXPRS = stratify.TRT,
                              OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = FALSE)
  g_CAN_regs_list[["gA"]] <- gA.sVars$regs

  gN.sVars <- process_regforms(regforms = gform.MONITOR, default.reg = gform.MONITOR.default, stratify.EXPRS = stratify.MONITOR,
                              OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = FALSE)
  g_CAN_regs_list[["gN"]] <- gN.sVars$regs

  # ------------------------------------------------------------------------------------------
  # DEFINE a single regression class
  # Perform S3 method dispatch on ALL_g_regs, which will determine the nested tree of SummaryModel objects
  # Perform fit and prediction
  # ------------------------------------------------------------------------------------------
  ALL_g_regs <- RegressionClass$new(RegressionForms = g_CAN_regs_list)
  ALL_g_regs$S3class <- "generic"
  modelfits.g0 <- newsummarymodel(reg = ALL_g_regs, DataStorageClass.g0 = OData)

  # modelfits.g0$fit(data = OData)
  modelfits.g0$fit(data = OData, predict = TRUE)

  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0

  # h_gN <- modelfits.g0$predictAeqa(newdata = OData)
  h_gN <- modelfits.g0$predictAeqa(n = OData$nobs)

  # ------------------------------------------------------------------------------------------
  # Evaluate indicator EVENT_IND that the person had experienced the outcome = 1 at any time of the follow-up:
  # ..... NOT REALLY NEEDED ......
  # ------------------------------------------------------------------------------------------
  # EVENT_IND <- "Delta"
  # if (EVENT_IND %in% names(OData$dat.sVar)) OData$dat.sVar[,(EVENT_IND):=NULL]
  # OData$dat.sVar[,(EVENT_IND):=as.integer(any(get(OUTCOME) %in% 1)), by = eval(ID)]
  # ------------------------------------------------------------------------------------------
  # Evaluate the indicator that this person was right-censored at some point in time:
  # ..... NOT REALLY NEEDED ......
  # ------------------------------------------------------------------------------------------
  # CENS_IND <- "AnyCensored"
  # if (CENS_IND %in% names(OData$dat.sVar)) OData$dat.sVar[,(CENS_IND):=NULL]
  # # noCENS.cat <- 0L; CENS <- c("C")
  # OData$dat.sVar[, (CENS_IND) := FALSE, by = eval(ID)]
  # for (Cvar in CENS) {
  #   OData$dat.sVar[, (CENS_IND) := get(CENS_IND) | any(!get(Cvar) %in% c(eval(noCENS.cat),NA)), by = eval(ID)]
  # }
  # OData$dat.sVar[, (CENS_IND) := as.integer(get(CENS_IND))]
  # ------------------------------------------------------------------------------------------
  # Evaluate the total duration of the follow-up for each observation
  # ..... NOT REALLY NEEDED ......
  # ------------------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------------------
  # Observed likelihood of (A,C,N) at each t, based on fitted object models in object modelfits.g0
  # ------------------------------------------------------------------------------------------
  # get back g_CAN_regs_list:
  OData$modelfits.g0 <- modelfits.g0

  ALL_g_regs <- modelfits.g0$reg
  g_CAN_regs_list <- ALL_g_regs$RegressionForms

  OData$modelfit.gC <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gC")]]
  OData$modelfit.gA <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gA")]]
  OData$modelfit.gN <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gN")]]

  g0.A <- OData$modelfit.gA$getcumprodAeqa()
  g0.C <- OData$modelfit.gC$getcumprodAeqa()
  g0.N <- OData$modelfit.gN$getcumprodAeqa()

  OData$dat.sVar[, c("g0.A", "g0.C", "g0.N", "g0.CAN") := list(g0.A, g0.C, g0.N, g0.A*g0.C*g0.N)]
  # newdat <- OData$dat.sVar[, list("g0.A" = g0.A, "g0.C" = g0.C, "g0.N" = g0.N, "g0.CAN" = g0.A*g0.C*g0.N)]

  return(OData)
}

# ---------------------------------------------------------------------------------------
# - BLOCK 3: evaluate weights based gstar.TRT, gstar.MONITOR and observed propensity scores g0, the input is modelfits.g0 and OData object
# ---------------------------------------------------------------------------------------
# Requires specification of probabilities for regimens of interest (either as rule followers or as counterfactual indicators)
# The output is person-specific data with evaluated weights, wts.DT, only observation-times with non-zero weight are kept
# Can be one regimen per single run of this block, which are then combined into a list of output datasets with lapply.
# Alternative is to allow input with several rules/regimens, which are automatically combined into a list of output datasets.
# ---------------------------------------------------------------------------------------
#' @export
#modelfits.g0,
get_weights <- function(OData, gstar.TRT = NULL, gstar.MONITOR = NULL) {
  nodes <- OData$nodes
  # OData$dat.sVar[, c("g0.CAN.compare") := list(h_gN)] # should be identical to g0.CAN
  # ------------------------------------------------------------------------------------------
  # Probabilities of counterfactual interventions under observed (A,C,N) at each t
  # Combine the propensity score for observed (g0.C, g0.A, g0.N) with the propensity scores for interventions (gstar.C, gstar.A, gstar.N):
  # ------------------------------------------------------------------------------------------------------------------------------
  # (1) gC.star: the indicator of not being censored.
  # (2) gA.star: prob of following one treatment rule; and
  # (3) gN.star prob following the monitoring regime; and
  # ------------------------------------------------------------------------------------------------------------------------------
  # indicator that the person is uncensored at each t (continuation of follow-up)
  # gstar.C <- "gstar.C"
  OData$dat.sVar[, "gstar.C" := as.integer(rowSums(.SD) == eval(OData$noCENS.cat)), .SDcols = nodes$Cnodes]
  # probability of following the rule at t, under intervention gstar.A on A(t)
  # **** NOTE ****
  # if gstar.TRT is a function then call it, if its a list of functions, then call one at a time.
  # if gstar.TRT returns more than one rule-column, estimate for each.
  if (!is.null(gstar.TRT)) {
    gstar.A <- as.name(gstar.TRT)
  } else {
    gstar.A <- as.name("g0.A") # use the actual observed exposure probability (no intervention on TRT)
  }
  # OData$dat.sVar[, "gstar.A" := get(gstar.A)]

  # probability of monitoring N(t)=1 under intervention gstar.N on N(t)
  # **** NOTE ****
  # if gstar.MONITOR is a function then call it, if its a list of functions, then call one at a time.
  # if gstar.MONITOR returns more than one rule-column, use each.
  if (!is.null(gstar.MONITOR)) {
    gstar.N <- as.name(gstar.MONITOR)
  } else {
    gstar.N <- as.name("g0.N") # use the actual observed monitoring probability (no intervention on MONITOR)
  }
  # OData$dat.sVar[, "gstar.N" := get(gstar.N)]

  # Joint probability for all 3:
  OData$dat.sVar[, "gstar.CAN" := gstar.C * eval(gstar.A) * eval(gstar.N)]

  # Weights by time and cummulative weights by time:
  OData$dat.sVar[, "wt.by.t" := gstar.CAN / g0.CAN, by = eval(nodes$IDnode)][, "cumm.IPAW" := cumprod(wt.by.t), by = eval(nodes$IDnode)]
  # OData$dat.sVar[1:100,]

  # -------------------------------------------------------------------------------------------
  # Weight stabilization - get emp P(followed rule at time t | followed rule up to now)
  # -------------------------------------------------------------------------------------------
  nIDs <- OData$nuniqueIDs
  # THE ENUMERATOR: the total sum of subjects followed the rule gstar.A at t
  # THE DENOMINATOR: divide above by the total number of subjects who were still at risk of NOT FOLLOWING the rule at t
  # i.e., followed rule at t-1, assume at the first time-point EVERYONE was following the rule (so denominator = n)
  # (The total sum of all subjects who WERE AT RISK at t)
  # In memory version (all at once inside data.table, by reference):
  # OData$dat.sVar[, N.follow.rule := sum(eval(gstar.A), na.rm = TRUE), by = eval(t.name)]
  # OData$dat.sVar[, cum.stab.P2 := cumprod(N.follow.rule / shift(N.follow.rule, fill = nIDs, type = "lag")), by = eval(ID)]

  # (MUCH FASTER) Version outside data.table, then merge back results:
  n.follow.rule.t <- OData$dat.sVar[, list(N.follow.rule = sum(eval(gstar.A), na.rm = TRUE)), by = eval(nodes$tnode)]
  n.follow.rule.t[, N.risk := shift(N.follow.rule, fill = nIDs, type = "lag")][, stab.P := N.follow.rule / N.risk][, cum.stab.P := cumprod(stab.P)]

  n.follow.rule.t[, c("N.risk", "stab.P") := list(NULL, NULL)]
  setkeyv(n.follow.rule.t, cols=nodes$tnode)
  OData$dat.sVar <- OData$dat.sVar[n.follow.rule.t, on = nodes$tnode]
  # equivalent: OData$dat.sVar <- merge(OData$dat.sVar, n.follow.rule.t, by = nodes$tnode)
  setkeyv(OData$dat.sVar, cols = c(nodes$IDnode, nodes$tnode))

  # Disabled: causes a bug (resuts in cumm.IPAW=0 when it shouldn't be)
  # OData$dat.sVar[cumm.IPAW < (10^-5), cum.stab.P := 0]
  # Disabled: remove all observation-times that got zero weight:
  # OData$dat.sVar[cumm.IPAW > 0, ]

  # multiply the weight by stabilization factor (numerator) (doesn't do anything for saturated MSMs, since they cancel):
  OData$dat.sVar[, cumm.IPAW := cum.stab.P * cumm.IPAW]

  Ynode <- nodes$Ynode # Get the outcome var:
  # shifted.OUTCOME <- Ynode%+%".tplus1"

  # OData$dat.sVar[, (shifted.OUTCOME) := shift(get(Ynode), n = 1L, type = "lead"), by = eval(nodes$IDnode)]
  # OData$dat.sVar <- OData$dat.sVar[!is.na(get(shifted.OUTCOME)), ] # drop and over-write previous data.table, removing last rows.

  # Multiply the shifted outcomes by the current (cummulative) weight cumm.IPAW:
  OData$dat.sVar[, "Wt.OUTCOME" := get(Ynode)*cumm.IPAW]
  # OData$dat.sVar[, "Wt.OUTCOME" := get(shifted.OUTCOME)*cumm.IPAW]
  # Row indices for all subjects at t who had the event at t+1 (NOT USING)
  # row_idx_outcome <- OData$dat.sVar[, .I[get(shifted.OUTCOME) %in% 1L], by = eval(ID)][["V1"]]
  # OData$dat.sVar[101:200, ]

  # Make a copy of the data.table only with relevant columns and keeping only the observations with non-zero weights
  # , na.rm = TRUE
  wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, "wt.by.t", "cumm.IPAW", "cum.stab.P", Ynode, "Wt.OUTCOME"), with = FALSE] # [wt.by.t > 0, ]
  # wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, "wt.by.t", "cumm.IPAW", "cum.stab.P", shifted.OUTCOME, "Wt.OUTCOME"), with = FALSE] # [wt.by.t > 0, ]
  wts.DT[, "rule.name.TRT" := eval(as.character(gstar.A))]
  wts.DT[, "rule.name.MONITOR" := eval(as.character(gstar.N))]

  # -------------------------------------------------------------------------------------------
  # NEED TO CLEAN UP OData$dat.sVar TO MAKE SURE ITS IN EXACTLY THE SAME STATE WHEN THIS FUNCTION WAS CALLED
  # -------------------------------------------------------------------------------------------
  # "g0.A", "g0.C", "g0.N", "g0.CAN",
  OData$dat.sVar[, c("gstar.C", "gstar.CAN", "Wt.OUTCOME", "cumm.IPAW", "wt.by.t", "cum.stab.P", "N.follow.rule") := list(NULL, NULL, NULL, NULL, NULL, NULL, NULL)]
  return(wts.DT)
}

# ---------------------------------------------------------------------------------------
# - BLOCK 4A: Non-parametric (saturated) MSM for survival, with weight stabilization,
# input either single weights dataset or a list of weights datasets.
# Each dataset containing weights non-zero weights for single regimen
# ---------------------------------------------------------------------------------------
#' @export
get_survNPMSM <- function(data.wts, OData) {
  nodes <- OData$nodes
  t.name <- nodes$tnode
  Ynode <- nodes$Ynode
  # shifted.OUTCOME <- as.name(Ynode%+%".tplus1")

  # CRUDE HAZARD ESTIMATE AND KM SURVIVAL:
  ht.crude <- data.wts[cumm.IPAW > 0, .(ht.KM = sum(eval(Ynode), na.rm = TRUE) / .N), by = eval(t.name)][, St.KM := cumprod(1 - ht.KM)]
  setkeyv(ht.crude, cols = t.name)

  # THE ENUMERATOR FOR THE HAZARD AT t: the weighted sum of subjects who had experienced the event at t:
  sum_Ywt <- data.wts[, .(sum_Y_IPAW = sum(Wt.OUTCOME, na.rm = TRUE)), by = eval(t.name)]; setkeyv(sum_Ywt, cols = t.name)
  # sum_Ywt <- OData$dat.sVar[, .(sum_Y_IPAW=sum(Wt.OUTCOME)), by = eval(t.name)]; setkeyv(sum_Ywt, cols=t.name)
  # THE DENOMINATOR FOR THE HAZARD AT t: The weighted sum of all subjects who WERE AT RISK at t:
  # (equivalent to summing cummulative weights cumm.IPAW by t)
  sum_Allwt <- data.wts[, .(sum_all_IPAW = sum(cumm.IPAW, na.rm = TRUE)), by = eval(t.name)]; setkeyv(sum_Allwt, cols = t.name)
  # sum_Allwt <- OData$dat.sVar[, .(sum_all_IPAW=sum(cumm.IPAW)), by = eval(t.name)]; setkeyv(sum_Allwt, cols=t.name)
  # EVALUATE THE DISCRETE HAZARD ht AND SURVIVAL St OVER t
  St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("St.IPTW") := .(cumprod(1 - ht))]
  # St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("m1ht", "St") := .(1-ht, cumprod(1-ht))]

  St_ht_IPAW <- merge(St_ht_IPAW, ht.crude, all=TRUE)
  return(list(IPW_estimates = data.frame(St_ht_IPAW)))
}


# ----------------------------------------------------------------------
# TO DO: 0. START A NEW CLASS THAT INHERITS FROM BinaryOutcomeModel SPECIFICALLY FOR RUNNING H20 MODELS
# TO DO: 1. MAKE THIS INTO A CALL TO BinaryOutcomeModel OR ITS CHILD
# TO DO: 2. SPLIT get_survMSM INTO ESTIMATION AND INFERENCE PARTS
# ----------------------------------------------------------------------
runglmMSM <- function(OData, wts.all.rules, all_dummies, Ynode, verbose) {
  # Generic prediction fun for logistic regression coefs, predicts P(A = 1 | X_mat)
  # Does not handle cases with deterministic Anodes in the original data.
  logispredict = function(m.fit, X_mat) {
    eta <- X_mat[,!is.na(m.fit$coef), drop = FALSE] %*% m.fit$coef[!is.na(m.fit$coef)]
    pAout <- match.fun(FUN = m.fit$linkfun)(eta)
    return(pAout)
  }

  if (getopt("fit.package") %in% c("h2o", "h2oglm")) {
    if (verbose) message("...fitting hazard MSM with h2o::h2o.glm...")
    loadframe_t <- system.time(
      MSM.designmat.H2O <- OData$fast.load.to.H2O(wts.all.rules,
                                                  saveH2O = FALSE,
                                                  destination_frame = "MSM.designmat.H2O")
    )
    if (verbose) { print("time to load the design mat into H2OFRAME: "); print(loadframe_t) }
    # OData$fast.load.to.H2O(wts.all.rules, )
    # temp.csv.path <- file.path(tempdir(), "wts.all.rules.csv~")
    # data.table::fwrite(wts.all.rules, temp.csv.path, turbo = TRUE)
    # MSM.designmat.H2O <- h2o::h2o.uploadFile(path = temp.csv.path, parse_type = "CSV", destination_frame = "MSM.designmat.H2O")

    m.fit_h2o <- try(h2o::h2o.glm(y = Ynode,
                                  x = all_dummies,
                                  intercept = FALSE,
                                  weights_column = "cumm.IPAW",
                                  training_frame = MSM.designmat.H2O,
                                  family = "binomial",
                                  standardize = FALSE,
                                  solver = c("IRLSM"), # solver = c("L_BFGS"),
                                  max_iterations = 50,
                                  lambda = 0L),
                silent = TRUE)

    out_coef <- vector(mode = "numeric", length = length(all_dummies))
    out_coef[] <- NA
    names(out_coef) <- c(all_dummies)
    out_coef[names(m.fit_h2o@model$coefficients)[-1]] <- m.fit_h2o@model$coefficients[-1]
    m.fit <- list(coef = out_coef, linkfun = "logit_linkinv", fitfunname = "h20")
    wts.all.rules[, glm.IPAW.predictP1 := as.vector(h2o::h2o.predict(m.fit_h2o, newdata = MSM.designmat.H2O)[,"p1"])]
  } else {
    if (verbose) message("...fitting hazard MSM with speedglm::speedglm.wfit...")
    Xdesign.mat <- as.matrix(wts.all.rules[, all_dummies, with = FALSE])
    m.fit <- try(speedglm::speedglm.wfit(
                                       X = Xdesign.mat,
                                       y = as.numeric(wts.all.rules[[Ynode]]),
                                       intercept = FALSE,
                                       family = binomial(),
                                       weights = wts.all.rules[["cumm.IPAW"]],
                                       trace = FALSE),
                        silent = TRUE)
    if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
      if (verbose) message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", m.fit)
      ctrl <- glm.control(trace = FALSE)
      SuppressGivenWarnings({
        m.fit <- stats::glm.fit(x = Xdesign.mat,
                                y = as.numeric(wts.all.rules[[Ynode]]),
                                family = binomial(),
                                intercept = FALSE, control = ctrl)
      }, GetWarningsToSuppress())
    }
    m.fit <- list(coef = m.fit$coef, linkfun = "logit_linkinv", fitfunname = "speedglm")
    if (verbose) {
      print("MSM fits"); print(m.fit$coef)
    }
    wts.all.rules[, glm.IPAW.predictP1 := logispredict(m.fit, Xdesign.mat)]
  }
  return(list(wts.all.rules = wts.all.rules, m.fit = m.fit))
}


# ---------------------------------------------------------------------------------------
# - BLOCK 4C: Parametric MSM pooling many regimens, includes weight stabilization and parametric MSM (can include saturated MSM)
# Alternative for MSM, takes a list of data sets with weights and MSMregform
# This block runs glm.fit got obtain MSM estimates (weighted regression) that uses the outcomes from many regimens, with dummy indicators for each input regime
# Might choose between two types of MSMs (sat vs. parametric) with S3 dispatch or by missing arg
# ---------------------------------------------------------------------------------------
#' @export
get_survMSM <- function(OData, data.wts.list, tjmin, tjmax, use.weights = TRUE, trunc.weights = Inf,
                        est.name = "IPAW", t.periods.RDs, verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t.name <- nodes$tnode
  Ynode <- nodes$Ynode
  # shifted.OUTCOME <- Ynode%+%".tplus1"

  # 2a. Stack the weighted data sets:
  wts.all.rules <- rbindlist(data.wts.list)
  rules.TRT <- sort(unique(wts.all.rules[["rule.name.TRT"]]))

  if (verbose) print("performing estimation for rules: " %+% paste(rules.TRT, collapse=","))

  # 2b. Remove all observations with 0 weights and run speedglm on design matrix with no intercept
  wts.all.rules <- wts.all.rules[!is.na(cumm.IPAW) & !is.na(eval(as.name(Ynode))) & (cumm.IPAW > 0), ]
  # 2c. If trunc.weights < Inf, do truncation of the weights
  if (trunc.weights < Inf) {
    wts.all.rules[cumm.IPAW > trunc.weights, cumm.IPAW := trunc.weights]
  }
  # 2d. use.weights==FALSE, do a crude estimator by setting all non-zero weights to 1
  if (use.weights == FALSE) {
    wts.all.rules[cumm.IPAW > 0, cumm.IPAW := 1L]
  }

  # 2.e. define all observed sequence of periods (t's)
  mint <- min(wts.all.rules[[t.name]])
  maxt <- max(wts.all.rules[[t.name]])
  periods <- mint:maxt
  periods_idx <- seq_along(periods)
  if (verbose) {
    print("periods"); print(periods)
  }

  # 3. Create the dummies I(d == gstar.TRT) for the logistic MSM for d-specific hazard
  all.d.dummies <- NULL
  for( dummy.j in rules.TRT ){
    wts.all.rules[, (dummy.j) := as.integer(rule.name.TRT %in% dummy.j)]
    all.d.dummies <- c(all.d.dummies, dummy.j)
  }
  # 4. Create the dummies I(t in interval.j), where interval.j defined by intervals of time of increasing length
  all.t.dummies <- NULL
  for( year.j in 1:length(tjmin)){
    dummy.j <- paste("Periods.",tjmin[year.j],"to",tjmax[year.j],sep="")
    wts.all.rules[, (dummy.j) := as.integer(eval(as.name(t.name)) >= tjmin[year.j] & eval(as.name(t.name)) <= tjmax[year.j])]
    all.t.dummies <- c(all.t.dummies, dummy.j)
  }
  # 5. Create interaction dummies I(t in interval.j & d == gstar.TRT)
  for (d.dummy in all.d.dummies) {
    for (t.dummy in all.t.dummies) {
      if (verbose) print(t.dummy %+% "_" %+% d.dummy)
      wts.all.rules[, (t.dummy %+% "_" %+% d.dummy) := as.integer(eval(as.name(t.dummy)) & eval(as.name(d.dummy)))]
    }
  }
  setkeyv(wts.all.rules, cols = c(nodes$IDnode, nodes$tnode))
  all_dummies <-  paste(sapply(all.d.dummies, function(x) {
                        return(paste(paste(paste(all.t.dummies, x, sep="_"), sep="")))
                        }))

  # 6. fit the hazard MSM
  resglmMSM <- runglmMSM(OData, wts.all.rules, all_dummies, Ynode, verbose)
  wts.all.rules <- resglmMSM$wts.all.rules
  m.fit <- resglmMSM$m.fit

  #### For variable estimation, GET IC and SE FOR BETA's
  beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                            MSMdata = wts.all.rules, MSMdesign = as.matrix(wts.all.rules[, all_dummies, with = FALSE]),
                            MSMpredict = "glm.IPAW.predictP1", IPW_MSMestimator = use.weights)

  # 7. Compute the Survival curves under each d
  S2.IPAW <- hazard.IPAW <- rep(list(rep(NA,maxt-mint+1)), length(rules.TRT))
  design.t.d <- rep(list(matrix(0L, ncol = length(all_dummies), nrow = length(mint:maxt))), length(rules.TRT))
  IC.Var.S.d <- vector(mode = "list", length(rules.TRT))
  names(S2.IPAW) <- names(hazard.IPAW) <- names(design.t.d) <- names(IC.Var.S.d) <- rules.TRT

  # the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
  # the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
  # p.coef - number of time-specific coefficients in the MSM
  p.coef <- length(tjmin)
  # periods <- mint : maxt
  design.t <- matrix(0L, ncol = p.coef, nrow = length(periods))
  for (t.indx in seq_along(periods)) {
    col.idx <- which(tjmin <= periods[t.indx] & tjmax >= periods[t.indx])
    design.t[t.indx, col.idx] <- 1
  }

  if (verbose) message("...evaluating survival & SEs based on MSM hazard fit and the estimated IC...")

  for(d.j in names(S2.IPAW)) {
    for(period.idx in seq_along(periods)){
      period.j <- periods[period.idx]
      rev.term <- paste0("Periods.",tjmin[max(which(tjmin<=period.j))],"to",tjmax[min(which(tjmax>=period.j))],"_",d.j)
      hazard.IPAW[[d.j]][period.idx] <- 1 / (1 + exp(-m.fit$coef[rev.term]))
      S2.IPAW[[d.j]][period.idx] <- (1-1/(1 + exp(-m.fit$coef[rev.term])))
    }

    S2.IPAW[[d.j]] <- cumprod(S2.IPAW[[d.j]])

    d.idx <- which(names(S2.IPAW) %in% d.j)
    set_cols <- seq((d.idx - 1) * ncol(design.t) + 1, (d.idx) * ncol(design.t))
    design.t.d[[d.j]][,set_cols] <- design.t

    #### GET IC and SE FOR Sd(t)
    # S.d.t.predict - MSM survival estimates for one regimen
    # h.d.t.predict - MSM hazard estimates for one regimen
    # design.d.t - d-specific matrix of dummy indicators for each t, i.e., d(m(t,d))/t
    # IC.O - observation-sepcific IC estimates for MSM coefs
    IC.Var.S.d[[d.j]] <- getSE.S(nID = nID,
                                 S.d.t.predict = S2.IPAW[[d.j]],
                                 h.d.t.predict = hazard.IPAW[[d.j]],
                                 design.d.t = design.t.d[[d.j]],
                                 IC.O = beta.IC.O.SEs[["IC.O"]])
  }

  output.MSM <- round(m.fit$coef,2)
  output.MSM <- cbind("Terms" = names(m.fit$coef), output.MSM)
  colnames(output.MSM) <- c("Terms",ifelse(trunc.weights == Inf && use.weights, "IPAW", ifelse(trunc.weights < Inf && use.weights, "truncated IPAW", "no weights")))
  rownames(output.MSM) <- NULL

  (quant99 <- quantile(wts.all.rules[["cumm.IPAW"]],p=0.99))
  (quant999 <- quantile(wts.all.rules[["cumm.IPAW"]],p=0.999))
  cutoffs <- c(0,0.5,1,10,20,30,40,50,100,150)
  IPAWdist <- makeSumFreqTable(table(wts.all.rules[["cumm.IPAW"]]),c(0,0.5,1,10,20,30,40,50,100,150),"Stabilized IPAW")

  ## RD:
  getSE_table_d_by_d <- function(S2.IPAW, IC.Var.S.d, nID, t.period.val.idx) {
    se.RDscale.Sdt.K <- matrix(NA, nrow = length(S2.IPAW), ncol = length(S2.IPAW))
    colnames(se.RDscale.Sdt.K) <- names(S2.IPAW)
    rownames(se.RDscale.Sdt.K) <- names(S2.IPAW)
    for (d1.idx in seq_along(names(S2.IPAW))) {
      for (d2.idx in seq_along(names(S2.IPAW))) {
        #### GET SE FOR RD(t)=Sd1(t) - Sd2(t)
        se.RDscale.Sdt.K[d1.idx, d2.idx] <- getSE.RD.d1.minus.d2(nID = nID,
                                                                IC.S.d1 = IC.Var.S.d[[d1.idx]][["IC.S"]],
                                                                IC.S.d2 = IC.Var.S.d[[d2.idx]][["IC.S"]])[t.period.val.idx]
      }
    }
    return(se.RDscale.Sdt.K)
  }

  RDs.IPAW.tperiods <- vector(mode = "list", length = length(t.periods.RDs))
  names(RDs.IPAW.tperiods) <- "RDs_for_t" %+% t.periods.RDs
  for (t.idx in seq(t.periods.RDs)) {
    # t.period.val <- t.periods.RDs[t.period.val.idx]
    t.period.val.idx <- periods_idx[periods %in% t.periods.RDs[t.idx]]
    se.RDscale.Sdt.K <- getSE_table_d_by_d(S2.IPAW, IC.Var.S.d, nID, t.period.val.idx)
    RDs.IPAW.tperiods[[t.idx]] <- make.table.m0(S2.IPAW, RDscale = TRUE, t.period = t.period.val.idx, nobs = nrow(wts.all.rules), esti = est.name, se.RDscale.Sdt.K = se.RDscale.Sdt.K)
  }

  return(list(St = S2.IPAW, MSM.fit = m.fit, output.MSM = output.MSM,
              IPAWdist = IPAWdist, RDs.IPAW.tperiods = RDs.IPAW.tperiods))
}
