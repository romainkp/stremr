
# ------------------------------------------------------------------
# - BLOCK 1: Process inputs and define OData R6 object
# ------------------------------------------------------------------
#' @export
get_Odata <- function(data, ID = "Subj_ID", t.name = "time_period", covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                      noCENS.cat = 0L, verbose = getOption("stremr.verbose")) {
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
  modelfits.g0$fit(data = OData, predict = TRUE)
  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0
  h_gN <- modelfits.g0$predictAeqa(n = OData$nobs)

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
  OData$dat.sVar[, "gstar.C" := as.integer(rowSums(.SD, na.rm = TRUE) == eval(OData$noCENS.cat)), .SDcols = nodes$Cnodes]
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

  # Joint probability for all 3:
  OData$dat.sVar[, "gstar.CAN" := gstar.C * eval(gstar.A) * eval(gstar.N)]
  # Weights by time and cummulative weights by time:
  OData$dat.sVar[, "wt.by.t" := gstar.CAN / g0.CAN, by = eval(nodes$IDnode)][, "cumm.IPAW" := cumprod(wt.by.t), by = eval(nodes$IDnode)]

  # -------------------------------------------------------------------------------------------
  # Weight stabilization - get emp P(followed rule at time t | followed rule up to now)
  # -------------------------------------------------------------------------------------------
  nIDs <- OData$nuniqueIDs
  # THE ENUMERATOR: the total sum of subjects followed the rule gstar.A at t
  # THE DENOMINATOR: divide above by the total number of subjects who were still at risk of NOT FOLLOWING the rule at t
  # i.e., followed rule at t-1, assume at the first time-point EVERYONE was following the rule (so denominator = n)
  # (The total sum of all subjects who WERE AT RISK at t)
  # (FASTER) Version outside data.table, then merge back results:
  n.follow.rule.t <- OData$dat.sVar[, list(N.follow.rule = sum(eval(gstar.A), na.rm = TRUE)), by = eval(nodes$tnode)]
  n.follow.rule.t[, N.risk := shift(N.follow.rule, fill = nIDs, type = "lag")][, stab.P := ifelse(N.risk > 0, N.follow.rule / N.risk, 0)][, cum.stab.P := cumprod(stab.P)]
  n.follow.rule.t[, c("N.risk", "stab.P") := list(NULL, NULL)]
  setkeyv(n.follow.rule.t, cols = nodes$tnode)
  OData$dat.sVar <- OData$dat.sVar[n.follow.rule.t, on = nodes$tnode]
  setkeyv(OData$dat.sVar, cols = c(nodes$IDnode, nodes$tnode))
  # Disabled: remove all observation-times that got zero weight:
  # OData$dat.sVar <- OData$dat.sVar[cumm.IPAW > 0, ]

  # multiply the weight by stabilization factor (numerator) (doesn't do anything for saturated MSMs, since they cancel):
  OData$dat.sVar[, cumm.IPAW := cum.stab.P * cumm.IPAW]
  Ynode <- nodes$Ynode # Get the outcome var:
  # Multiply the shifted outcomes by the current (cummulative) weight cumm.IPAW:
  OData$dat.sVar[, "Wt.OUTCOME" := get(Ynode)*cumm.IPAW]
  # Row indices for all subjects at t who had the event at t+1 (NOT USING)
  # row_idx_outcome <- OData$dat.sVar[, .I[get(Ynode) %in% 1L], by = eval(ID)][["V1"]]

  # Make a copy of the data.table only with relevant columns and keeping only the observations with non-zero weights
  wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, "wt.by.t", "cumm.IPAW", "cum.stab.P", Ynode, "Wt.OUTCOME"), with = FALSE] # [wt.by.t > 0, ]
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
get_survNPMSM <- function(wts.data, OData) {
  nodes <- OData$nodes
  t.name <- nodes$tnode
  Ynode <- nodes$Ynode
  # CRUDE HAZARD ESTIMATE AND KM SURVIVAL:
  ht.crude <- wts.data[cumm.IPAW > 0, .(ht.KM = sum(eval(as.name(Ynode)), na.rm = TRUE) / .N), by = eval(t.name)][, St.KM := cumprod(1 - ht.KM)]
  setkeyv(ht.crude, cols = t.name)
  # THE ENUMERATOR FOR THE HAZARD AT t: the weighted sum of subjects who had experienced the event at t:
  sum_Ywt <- wts.data[, .(sum_Y_IPAW = sum(Wt.OUTCOME, na.rm = TRUE)), by = eval(t.name)]; setkeyv(sum_Ywt, cols = t.name)
  # sum_Ywt <- OData$dat.sVar[, .(sum_Y_IPAW=sum(Wt.OUTCOME)), by = eval(t.name)]; setkeyv(sum_Ywt, cols=t.name)
  # THE DENOMINATOR FOR THE HAZARD AT t: The weighted sum of all subjects who WERE AT RISK at t:
  # (equivalent to summing cummulative weights cumm.IPAW by t)
  sum_Allwt <- wts.data[, .(sum_all_IPAW = sum(cumm.IPAW, na.rm = TRUE)), by = eval(t.name)]; setkeyv(sum_Allwt, cols = t.name)
  # sum_Allwt <- OData$dat.sVar[, .(sum_all_IPAW=sum(cumm.IPAW)), by = eval(t.name)]; setkeyv(sum_Allwt, cols=t.name)
  # EVALUATE THE DISCRETE HAZARD ht AND SURVIVAL St OVER t
  St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("St.IPTW") := .(cumprod(1 - ht))]
  # St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("m1ht", "St") := .(1-ht, cumprod(1-ht))]
  St_ht_IPAW <- merge(St_ht_IPAW, ht.crude, all=TRUE)
  return(list(IPW_estimates = data.frame(St_ht_IPAW)))
}

# ---------------------------------------------------------------------------------------
#' Estimate Survival with a particular MSM for the survival-hazard function using previously fitted weights.
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the treatment/exposure and monitoring processes based on
#' the user-specified Marginal Structural Model (MSM) for the counterfactual survival function.
#'
#' This routine will run the weighted logistic regression using the (possibly-weighted) outcomes from many regimens, with dummy indicators for each treatment/monitoring
#' regimen available in \code{wts.data} and each follow-up time interval specified in \code{t.breaks}.
#' When \code{use.weights = TRUE}, the logistic regression for the survival hazard is weighted by the \strong{IPW} (Inverse Probability-Weighted or Horvitz-Thompson) estimated weights
#' in \code{wts.data}. These IPW weights are based on the previously fitted propensity scores (function \code{get_fits}), allowing
#' adjustment for confounding by possibly non-random assignment to exposure and monitoring and possibly informative right-censoring.
#' @param OData The object returned by function \code{get_fits}. Contains the input dat and the previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param wts.data A list of \code{data.table}s, each data set is a result of calling the function \code{get_weights}. Must contain the treatment/monitoring rule-specific estimated IPTW weights.
#' This argument can be also a single \code{data.table} obtained with \code{data.table::rbindlist(wts.data)}.
#' @param t.breaks The vector of integer (or numeric) breaks that defines the dummy indicators of the follow-up in the observed data. Used for fitting the parametric (or saturated) MSM for the survival hazard. See Details.
#' @param use.weights Logical value. Set to \code{FALSE} to ignore the weights in \code{wts.data} and fit a "crude" MSM that does not adjust for the possible confounding due to non-random
#' assignment of the exposure/censoring and monitoring indicators.
#' @param trunc.weights Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc.weights} will be truncated.
#' @param est.name A string naming the current MSM estimator. Ignored by the current routine but is used when generating reports with \code{make_report_rmd}.
#' @param getSEs A logical indicator. Set to \code{TRUE} to evaluate the standard errors for the estimated survival by using the MSM influence curve.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#'
#' @section Details:
#' **********************************************************************
#' \code{t.breaks} is used for defining the time-intervals of the MSM coefficients for estimation of the survival hazard function.
#' The first value in \code{t.breaks} defines a dummy variable (indicator) for a fully closed interval, with each subsquent value in \code{t.breaks} defining a single right-closed time-interval.
#' For example, \code{t.breaks = c(0,1)} will define the MSM dummy indicators: I(min(t) <= t <=0 ), I(0 < t <= 1) and I(1 < t <= max(t)).
#' On the other hand \code{t.breaks = c(1)} will define the following (more parametric) MSM dummy indicators: I(min(t) <= t <=1 ) and I(1 < t <= max(t)).
#' If omitted, the default is to define a saturated (non-parametric) MSM with a separate dummy variable for every unique period in the observed data.

#' @section MSM for the hazard:
#' **********************************************************************
#'
#' @return A named list with items containing the MSM estimation results:
#'  \itemize{
#'  \item \code{est.name} - .
#'  \item \code{St} - .
#'  \item \code{ht} - .
#'  \item \code{MSM.fit} - .
#'  \item \code{MSM.intervals} - .
#'  \item \code{IC.Var.S.d} - .
#'  \item \code{nID} - .
#'  \item \code{wts.data} - .
#'  \item \code{use.weights} - .
#'  \item \code{trunc.weights} - .
#' }
#' @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/4_get_survMSM_example.R
#' @export
get_survMSM <- function(OData, wts.data, t.breaks, use.weights = TRUE, trunc.weights = 10^6, est.name = "IPAW", getSEs = TRUE, verbose = getOption("stremr.verbose")) {
  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t.name <- nodes$tnode
  Ynode <- nodes$Ynode

  # 2a. Stack the weighted data sets, if those came in a list:
  if (is.data.table(wts.data)) {
    # ...do nothing...
  } else if (is.list(wts.data)) {
    assert_that(all(sapply(wts.data, is.data.table)))
    wts.data <- rbindlist(wts.data)
  } else {
    stop("...wts.data arg must be either a list of rule-specific weight data.tables or one combined weight data.table...")
  }

  rules.TRT <- sort(unique(wts.data[["rule.name.TRT"]]))
  if (verbose) print("performing estimation for found TRT rules in column 'rule.name.TRT': " %+% paste(rules.TRT, collapse=","))

  # 2b. Remove all observations with 0 weights and run speedglm on design matrix with no intercept
  wts.data <- wts.data[!is.na(cumm.IPAW) & !is.na(eval(as.name(Ynode))) & (cumm.IPAW > 0), ]
  setkeyv(wts.data, cols = c(nodes$IDnode, nodes$tnode))

  # 2c. If trunc.weights < Inf, do truncation of the weights
  if (trunc.weights < Inf) {
    wts.data[cumm.IPAW > trunc.weights, cumm.IPAW := trunc.weights]
  }
  # 2d. use.weights==FALSE, do a crude estimator by setting all non-zero weights to 1
  if (use.weights == FALSE) {
    wts.data[cumm.IPAW > 0, cumm.IPAW := 1L]
  }

  # 2.e. define all observed sequence of periods (t's)
  mint <- min(wts.data[[t.name]], na.rm = TRUE); maxt <- max(wts.data[[t.name]], na.rm = TRUE)
  periods <- (mint:maxt)
  periods_idx <- seq_along(periods)
  if (verbose) {
    print("periods"); print(periods)
  }

  # 2.f. default t.breaks, error checks for t.breaks, plus padding w/ mint & maxt:
  if (missing(t.breaks)) {
    # default t.breaks is to use a saturated (non-parametric) MSM
    t.breaks <- sort(periods)
    if (verbose)
      message("running with default 't.breaks': (" %+%
        paste0(t.breaks, collapse = ",") %+%
        "); \nNote: such 't.breaks' define a separate coefficient for every unique follow-up time period resulting in a saturated (non-parametric) MSM.")
  }
  if (length(unique(t.breaks)) < length(t.breaks)) {
    stop("all t.breaks must be unique")
  }
  if (!all(t.breaks %in% periods)) {
    stop("all t.breaks must be contained between minimum and maximum follow-up periods:" %+% t.breaks[!(t.breaks %in% periods)])
  }

  if (max(t.breaks) < maxt) t.breaks <- sort(c(t.breaks, maxt)) # pad on the right (if needed with maxt):

  # 3. Create the dummies I(d == gstar.TRT) for the logistic MSM for d-specific hazard
  all.d.dummies <- NULL
  for( dummy.j in rules.TRT ){
    wts.data[, (dummy.j) := as.integer(rule.name.TRT %in% dummy.j)]
    all.d.dummies <- c(all.d.dummies, dummy.j)
  }

  # 4. Create the dummies I(t in interval.j), where interval.j defined by intervals of time of increasing length
  all.t.dummies <- NULL
  t.breaks.mint <- c(mint, t.breaks) # pad t.breaks on the left (with mint)
  MSM.intervals <- matrix(NA, ncol = 2, nrow = length(t.breaks)) # save the actual intervals
  colnames(MSM.intervals) <- c("min.t", "max.t")
  t.per.inteval <- vector(mode = "list", length = nrow(MSM.intervals)) # save the vector of period vals that belong to each interval
  for (t.idx in 2:length(t.breaks.mint)) {
    low.t <- t.breaks.mint[t.idx - 1]
    high.t <- t.breaks.mint[t.idx]
    # First interval needs to be closed on both sides (includes the smallest obesrved follow-up, mint)
    if (t.idx == 2L) {
      dummy.j <- paste("Periods.", low.t, "to", high.t, sep="")
      MSM.intervals[t.idx - 1, ] <- c(low.t, high.t); t.per.inteval[[t.idx - 1]] <- unique(low.t:high.t)
      wts.data[, (dummy.j) := as.integer(eval(as.name(t.name)) >= low.t & eval(as.name(t.name)) <= high.t)]
    } else {
      dummy.j <- paste("Periods.", (low.t + 1), "to", high.t, sep="")
      MSM.intervals[t.idx - 1, ] <- c(low.t + 1, high.t); t.per.inteval[[t.idx - 1]] <- unique((low.t+1):high.t)
      wts.data[, (dummy.j) := as.integer(eval(as.name(t.name)) >= (low.t + 1) & eval(as.name(t.name)) <= high.t)]
    }
    print("defined t.dummy: " %+% dummy.j)
    all.t.dummies <- c(all.t.dummies, dummy.j)
  }

  # 5. Create interaction dummies I(t in interval.j & d == gstar.TRT)
  for (d.dummy in all.d.dummies) {
    for (t.dummy in all.t.dummies) {
      if (verbose)
        print(t.dummy %+% "_" %+% d.dummy)
      wts.data[, (t.dummy %+% "_" %+% d.dummy) := as.integer(eval(as.name(t.dummy)) & eval(as.name(d.dummy)))]
    }
  }

  all_dummies <-  paste(sapply(all.d.dummies, function(x) {
                        return(paste(paste(paste(all.t.dummies, x, sep="_"), sep="")))
                        }))

  # 6. fit the hazard MSM
  resglmMSM <- runglmMSM(OData, wts.data, all_dummies, Ynode, verbose)
  wts.data <- resglmMSM$wts.data
  m.fit <- resglmMSM$m.fit

  # 7. Compute the Survival curves under each d
  S2.IPAW <- hazard.IPAW <- rep(list(rep(NA,maxt-mint+1)), length(rules.TRT))
  names(S2.IPAW) <- names(hazard.IPAW) <- rules.TRT

  if (verbose) message("...evaluating MSM-based survival curves...")
  for(d.j in names(S2.IPAW)) {
    for(period.idx in seq_along(periods)) {
      period.j <- periods[period.idx] # the period of the follow-up for which we want to evaluate the MSM-based survival:
      # the dummy coefficient of the MSM that includes this time-point (period)
      # that is, find the appropriate right-closed interval from MSM.intervals matrix for a given period.j:
      int_idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
      if (!(period.j %in% t.per.inteval[[int_idx]])) stop("error while finding the appropriate MSM coefficient")
      d.j.idx <- which(all.d.dummies %in% d.j)
      MSM.term <- all_dummies[length(all.t.dummies)*(d.j.idx - 1) + int_idx]
      print("fetching MSM coefficient for period " %+% period.j %+% " and rule " %+% d.j %+% ": " %+% MSM.term)
      hazard.IPAW[[d.j]][period.idx] <- 1 / (1 + exp(-m.fit$coef[MSM.term]))
      S2.IPAW[[d.j]][period.idx] <- (1-1/(1 + exp(-m.fit$coef[MSM.term])))
    }
    S2.IPAW[[d.j]] <- cumprod(S2.IPAW[[d.j]])
  }

  if (verbose) message("...evaluating SEs based on MSM hazard fit and the estimated IC...")
  #### For variance (SEs), GET IC and SE FOR BETA's
  #### GET IC and SE FOR Sd(t)
  # S.d.t.predict - MSM survival estimates for one regimen
  # h.d.t.predict - MSM hazard estimates for one regimen
  # design.d.t - d-specific matrix of dummy indicators for each t, i.e., d(m(t,d))/t
  # IC.O - observation-specific IC estimates for MSM coefs
  if (getSEs) {
    design.d.t <- rep(list(matrix(0L, ncol = length(all_dummies), nrow = length(mint:maxt))), length(rules.TRT))
    IC.Var.S.d <- vector(mode = "list", length(rules.TRT))
    names(design.d.t) <- names(IC.Var.S.d) <- rules.TRT
    # the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
    # the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
    # p.coef - number of time-specific coefficients in the MSM
    p.coef <- length(tjmin)
    design.t <- matrix(0L, ncol = p.coef, nrow = length(periods))
    for (period.idx in seq_along(periods)) {
      period.j <- periods[period.idx]
      col.idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
      design.t[period.idx, col.idx] <- 1
    }
    beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                              MSMdata = wts.data, MSMdesign = as.matrix(wts.data[, all_dummies, with = FALSE]),
                              MSMpredict = "glm.IPAW.predictP1", IPW_MSMestimator = use.weights)

    for(d.j in names(S2.IPAW)) {
      d.idx <- which(names(S2.IPAW) %in% d.j)
      set_cols <- seq((d.idx - 1) * ncol(design.t) + 1, (d.idx) * ncol(design.t))
      design.d.t[[d.j]][,set_cols] <- design.t
      IC.Var.S.d[[d.j]] <- getSE.S(nID = nID,
                                   S.d.t.predict = S2.IPAW[[d.j]],
                                   h.d.t.predict = hazard.IPAW[[d.j]],
                                   design.d.t = design.d.t[[d.j]],
                                   IC.O = beta.IC.O.SEs[["IC.O"]])
    }
  } else {
    IC.Var.S.d <- NULL
  }
  MSM_out <- list(
              est.name = est.name,
              St = S2.IPAW,
              ht = hazard.IPAW,
              MSM.fit = m.fit,
              MSM.intervals = MSM.intervals,
              IC.Var.S.d = IC.Var.S.d,
              nID = nID, periods = periods,
              wts.data = wts.data,
              use.weights = use.weights,
              trunc.weights = trunc.weights
            )

  return(MSM_out)
}

