
# ------------------------------------------------------------------
# - BLOCK 1: Process inputs and define OData R6 object
# ------------------------------------------------------------------
#' @export
get_Odata <- function(data, ID = "Subj_ID", t = "time_period", covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                      noCENS.cat = 0L, verbose = FALSE) {
  gvars$verbose <- TRUE
  gvars$noCENS.cat <- noCENS.cat
  # if (verbose) {
    message("Running with the following setting: ");
    str(gvars$opts)
  # }
  if (missing(covars)) { # define time-varing covars (L) as everything else in data besides these vars
    covars <- setdiff(colnames(data), c(ID, t, CENS, TRT, MONITOR, OUTCOME))
  }
  # The ordering of variables in this list is the assumed temporal order!
  nodes <- list(Lnodes = covars, Cnodes = CENS, Anodes = TRT, Nnodes = MONITOR, Ynode = OUTCOME, IDnode = ID, tnode = t)
  OData <- DataStorageClass$new(Odata = data, nodes = nodes, noCENS.cat = noCENS.cat)
  # --------------------------------------------------------------------------------------------------------
  # Create dummies for factors
  # --------------------------------------------------------------------------------------------------------
  factor.Ls <- unlist(lapply(OData$dat.sVar, is.factor))
  factor.Ls <- names(factor.Ls)[factor.Ls]
  new.factor.names <- vector(mode="list", length=length(factor.Ls))
  names(new.factor.names) <- factor.Ls
  if (length(factor.Ls)>0)
    message("found factors in the data, these are being converted to binary indicators (first level excluded): " %+% paste0(factor.Ls, collapse=","))
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
  # ---------------------------------------------------------------------------
  # DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
  # Might expand this in the future to allow defining arbitrary summaries
  # ---------------------------------------------------------------------------
  lagnodes <- c(nodes$Cnodes, nodes$Anodes, nodes$Nnodes)
  newVarnames <- lagnodes %+% ".tminus1"
  OData$dat.sVar[, (newVarnames) := shift(.SD, n=1L, fill=1L, type="lag"), by=get(nodes$ID), .SDcols=(lagnodes)]

  # -------------------------------------------------------------------------------------------
  # Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
  # -------------------------------------------------------------------------------------------
  Ynode <- nodes$Ynode
  shifted.OUTCOME <- Ynode%+%".tplus1"
  OData$dat.sVar[, (shifted.OUTCOME) := shift(get(Ynode), n = 1L, type = "lead"), by = eval(nodes$IDnode)]
  # OData$dat.sVar <- OData$dat.sVar[!is.na(get(shifted.OUTCOME)), ] # drop and over-write previous data.table, removing last rows.

  for (Cnode in nodes$Cnodes) CheckVarNameExists(OData$dat.sVar, Cnode)
  for (Anode in nodes$Anodes) CheckVarNameExists(OData$dat.sVar, Anode)
  for (Nnode in nodes$Nnodes) CheckVarNameExists(OData$dat.sVar, Nnode)
  for (Ynode in nodes$Ynode)  CheckVarNameExists(OData$dat.sVar, Ynode)
  for (Lnode in nodes$Lnodes) CheckVarNameExists(OData$dat.sVar, Lnode)

  OData.R6 <- OData
  return(OData.R6)
}

# ------------------------------------------------------------------
# - BLOCK 2: define regression models, define a single RegressionClass & fit the propensity score for observed data, summary.g0 g0 (C,A,N)
# ------------------------------------------------------------------
#' @export
get_fits <- function(OData, gform.CENS, gform.TRT, gform.MONITOR,
                    stratify.CENS = NULL, stratify.TRT = NULL, stratify.MONITOR = NULL) {
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
  modelfits.g0$fit(data = OData)
  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0
  h_gN <- modelfits.g0$predictAeqa(newdata = OData)
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
  # OData$modelfit.gC$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$getPsAsW.models()
  OData$modelfit.gA <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gA")]]
  # OData$modelfit.gA$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$getPsAsW.models()
  OData$modelfit.gN <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gN")]]
  # OData$modelfit.gN$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$getfit

  g0.A <- OData$modelfit.gA$getcumprodAeqa()
  g0.C <- OData$modelfit.gC$getcumprodAeqa()
  g0.N <- OData$modelfit.gN$getcumprodAeqa()

  OData$dat.sVar[, c("g0.A", "g0.C", "g0.N", "g0.CAN") := list(g0.A, g0.C, g0.N, g0.A*g0.C*g0.N)]
  # newdat <- OData$dat.sVar[, list("g0.A" = g0.A, "g0.C" = g0.C, "g0.N" = g0.N, "g0.CAN" = g0.A*g0.C*g0.N)]

  return(OData)
  # return(list(IPW_estimates = data.frame(St_ht_IPAW), dataDT = OData$dat.sVar, modelfits.g0.R6 = modelfits.g0, OData.R6 = OData))
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
  # OData$dat.sVar[, N.follow.rule := sum(eval(gstar.A), na.rm = TRUE), by = eval(t)]
  # OData$dat.sVar[, cum.stab.P2 := cumprod(N.follow.rule / shift(N.follow.rule, fill = nIDs, type = "lag")), by = eval(ID)]

  # (MUCH FASTER) Version outside data.table, then merge back results:
  n.follow.rule.t <- OData$dat.sVar[, list(N.follow.rule = sum(eval(gstar.A), na.rm = TRUE)), by = eval(nodes$tnode)]
  n.follow.rule.t[, N.risk := shift(N.follow.rule, fill = nIDs, type = "lag")][, stab.P := N.follow.rule / N.risk][, cum.stab.P := cumprod(stab.P)]
  n.follow.rule.t[, c("N.risk", "stab.P"):=list(NULL, NULL)]
  setkeyv(n.follow.rule.t, cols=nodes$tnode)
  OData$dat.sVar <- OData$dat.sVar[n.follow.rule.t, on = nodes$tnode]
  # equivalent: OData$dat.sVar <- merge(OData$dat.sVar, n.follow.rule.t, by = nodes$tnode)
  OData$dat.sVar[cumm.IPAW < (10^-5), cum.stab.P := 0]
  setkeyv(OData$dat.sVar, cols = c(nodes$IDnode, nodes$tnode))
  # remove all observation-times that got zero weight:
  # OData$dat.sVar[cumm.IPAW > 0, ]
  # multiply the weight by stabilization factor (numerator) (doesn't do anything, since it cancels):
  OData$dat.sVar[, cumm.IPAW := cum.stab.P * cumm.IPAW]

  # Drop all observations with NA outcome:
  Ynode <- nodes$Ynode
  shifted.OUTCOME <- Ynode%+%".tplus1"
  # OData$dat.sVar[, (shifted.OUTCOME) := shift(get(Ynode), n = 1L, type = "lead"), by = eval(nodes$IDnode)]
  # OData$dat.sVar <- OData$dat.sVar[!is.na(get(shifted.OUTCOME)), ] # drop and over-write previous data.table, removing last rows.

  # Multiply the shifted outcomes by the current (cummulative) weight cumm.IPAW:
  OData$dat.sVar[, "Wt.OUTCOME" := get(shifted.OUTCOME)*cumm.IPAW]
  # Row indices for all subjects at t who had the event at t+1 (NOT USING)
  # row_idx_outcome <- OData$dat.sVar[, .I[get(shifted.OUTCOME) %in% 1L], by = eval(ID)][["V1"]]
  # OData$dat.sVar[101:200, ]

  # Make a copy of the data.table only with relevant columns and keeping only the observations with non-zero weights
  # , na.rm = TRUE
  wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, "wt.by.t", "cumm.IPAW", "cum.stab.P", shifted.OUTCOME, "Wt.OUTCOME"), with = FALSE] # [wt.by.t > 0, ]
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
# - BLOCK 4A: Non-parametric MSM for survival, no weight stabilization, input either single weights dataset or a list of weights datasets,
# Each dataset containing weights non-zero weights for single regimen
# ---------------------------------------------------------------------------------------
#' @export
get_survNP <- function(data.wts, OData) {
  nodes <- OData$nodes
  t <- nodes$tnode
  Ynode <- nodes$Ynode
  shifted.OUTCOME <- as.name(Ynode%+%".tplus1")

  # CRUDE HAZARD ESTIMATE AND KM SURVIVAL:
  ht.crude <- data.wts[cumm.IPAW > 0, .(ht.KM = sum(eval(shifted.OUTCOME), na.rm = TRUE) / .N), by = eval(t)][, St.KM := cumprod(1 - ht.KM)]
  setkeyv(ht.crude, cols = t)

  # THE ENUMERATOR FOR THE HAZARD AT t: the weighted sum of subjects who had experienced the event at t:
  sum_Ywt <- data.wts[, .(sum_Y_IPAW = sum(Wt.OUTCOME, na.rm = TRUE)), by = eval(t)]; setkeyv(sum_Ywt, cols = t)
  # sum_Ywt <- OData$dat.sVar[, .(sum_Y_IPAW=sum(Wt.OUTCOME)), by = eval(t)]; setkeyv(sum_Ywt, cols=t)
  # THE DENOMINATOR FOR THE HAZARD AT t: The weighted sum of all subjects who WERE AT RISK at t:
  # (equivalent to summing cummulative weights cumm.IPAW by t)
  sum_Allwt <- data.wts[, .(sum_all_IPAW = sum(cumm.IPAW, na.rm = TRUE)), by = eval(t)]; setkeyv(sum_Allwt, cols = t)
  # sum_Allwt <- OData$dat.sVar[, .(sum_all_IPAW=sum(cumm.IPAW)), by = eval(t)]; setkeyv(sum_Allwt, cols=t)
  # EVALUATE THE DISCRETE HAZARD ht AND SURVIVAL St OVER t
  St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("St.IPTW") := .(cumprod(1 - ht))]
  # St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("m1ht", "St") := .(1-ht, cumprod(1-ht))]

  St_ht_IPAW <- merge(St_ht_IPAW, ht.crude, all=TRUE)
  return(list(IPW_estimates = data.frame(St_ht_IPAW)))
}

# ---------------------------------------------------------------------------------------
# - BLOCK 4B: Saturated MSM pooling many regimens, includes weight stabilization and using closed-form soluaton for the MSM (can only do saturated MSM)
# ---------------------------------------------------------------------------------------
# This block runs an closed-form MSM regression for saturated MSM (dummy indicators for every interaction of t and regimen)
get_survMSM_1 <- function(data.wts.list, t) {
  ##############################################################
  ## Manual weighted logistic regression saturated MSM
  ##############################################################
    glmterm.names <- paste(sapply(all.d.dummies,function(x){
      return(paste(paste("I(",paste(all.t.dummies2,x,sep=" * "),")",sep="")))
    }))
    #### Fit manually
    NPIPAW.h2 <- rep(NA,2*16)
    names(NPIPAW.h2) <- glmterm.names
    NPcrude.h2 <- NPIPAW.h2
    count <- 0
    for(dtheta in c("dlow","dhigh")) {
      for(t in 1:16){
        count <- count+1
        relevantObs <- (long.IPAW.data[,"t"]%in%t & long.IPAW.data[,dtheta]%in%1)
        NPIPAW.h2[count] <- sum(long.IPAW.data[relevantObs, Yname] * long.IPAW.data[relevantObs, "IPAW"]) / sum(long.IPAW.data[relevantObs, "IPAW"])
        NPcrude.h2[count] <- sum(long.IPAW.data[relevantObs, Yname]) / length(long.IPAW.data[relevantObs, "IPAW"])
      }
    }
    glm.IPAW.h.coef2 <- log(NPIPAW.h2 / (1 - NPIPAW.h2))
    glm.h.coef2 <- log(NPcrude.h2 / (1 - NPcrude.h2))

    ##############################################################
    ## Compute the Survival curves under each d
    ##############################################################
    (mint <- min(long.IPAW.data[,"t"]))
    (maxt <- max(long.IPAW.data[,"t"]))
    S2.IPAW <- S2 <- rep(list(rep(NA,maxt-mint+1)),2)
    names(S2)  <- names(S2.IPAW) <- names(S2.IPAWtrunc) <- c("dlow","dhigh")

    for(d.j in names(S2)) {
      for(period.j in mint:maxt) {
        rev.term <- paste("I(",paste("Periods.",tjmin2[max(which(tjmin2<=period.j))],"to",tjmax2[min(which(tjmax2>=period.j))],sep="")," * ",d.j,")",sep="")
        S2.IPAW[[d.j]][period.j] <- (1-1/(1+exp(-glm.IPAW.h.coef2[rev.term])))
        S2[[d.j]][period.j] <- (1-1/(1+exp(-glm.h.coef2[rev.term])))
      }
    }
    S2 <- lapply(S2,cumprod)
    S2.IPAW <- lapply(S2.IPAW,cumprod)
    nrow.long.IPAW.data <- nrow(long.IPAW.data)
    ## Need this to compute inference
    IC.data <- long.IPAW.data[ ,c("ID","t","theta",Yname,"IPAW",all.t.dummies2,rules)] # contains data for all obs compatible with at least one rule
    return(list("nID"=nID,"mint"=mint,"maxt"=maxt,"IPAWdist"=IPAWdist,"quant99"=quant99,"quant999"=quant999,"nabove20"=nabove20,"glm.IPAW.h.coef2"=glm.IPAW.h.coef2,"glm.IPAWtrunc.h.coef2"=glm.IPAWtrunc.h.coef2,"glm.h.coef2"=glm.h.coef2,"S2"=S2,"S2.IPAW"=S2.IPAW,"all.t.dummies2"=all.t.dummies2,"IC.data"=IC.data,"rules"=rules))
  # return(list(IPW_estimates = data.frame(St_ht_IPAW)))
}

# ---------------------------------------------------------------------------------------
# - BLOCK 4C: Parametric MSM pooling many regimens, includes weight stabilization and parametric MSM (can include saturated MSM)
# ---------------------------------------------------------------------------------------
# Alternative for MSM, takes a list of data sets with weights and MSMregform
# This block runs glm.fit got obtain MSM estimates (weighted regression) that uses the outcomes from many regimens, with dummy indicators for each input regime
# Might choose between two types of MSMs (sat vs. parametric) with S3 dispatch or by missing arg
get_survMSM_2 <- function(data.wts.list, t, MSMregform) {
  # NOTE IMPLEMENTED
  # ....
  # ....
  # return(list(IPW_estimates = data.frame(St_ht_IPAW)))
}
