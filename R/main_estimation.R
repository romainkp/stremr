# adding to appease CRAN check with non-standard eval in data.table:
utils::globalVariables(c("gstar.CAN", "g0.CAN", "wt.by.t", "rule.follower.gCAN", "new.TRT.gstar",
                          "N.risk", "N.follow.rule", "stab.P", "cum.stab.P", "cum.IPAW",
                          "rule.name", "glm.IPAW.predictP1", "St.KM", "Wt.OUTCOME", "ht", "ht.KM", "EIC_i_t0", "EIC_i_tplus"))

# ---------------------------------------------------------------------------------------
#' Import data, define various nodes, define dummies for factor columns and define OData R6 object
#'
#' @param data Input data in long format. Can be a \code{data.frame} or a \code{data.table} with named columns, containing the time-varying covariates (\code{covars}),
#'  the right-censoring event indicator(s) (\code{CENS}), the exposure variable(s) (\code{TRT}), the monitoring process variable(s) (\code{MONITOR})
#'  and the survival OUTCOME variable (\code{OUTCOME}).
#' @param ID Unique subject identifier column name in \code{data}.
#' @param t_name The name of the time/period variable in \code{data}.
#' @param covars Vector of names with time varying and baseline covariates in \code{data}. This argument does not need to be specified, by default all variables
#' that are not in \code{ID}, \code{t}, \code{CENS}, \code{TRT}, \code{MONITOR} and \code{OUTCOME} will be considered as covariates.
#' @param CENS Column name of the censoring variable(s) in \code{data}.
#' Each separate variable specified in \code{CENS} can be either binary (0/1 valued integer) or categorical (integer).
#' For binary indicators of CENSoring, the value of 1 indicates the CENSoring or end of follow-up event (this cannot be changed).
#' For categorical CENSoring variables, by default the value of 0 indicates no CENSoring / continuation of follow-up and other
#' values indicate different reasons for CENSoring.
#' Use the argument \code{noCENScat} to change the reference (continuation of follow-up) category from default 0 to any other value.
#' (NOTE: Changing \code{noCENScat} has zero effect on coding of the binary CENSoring variables, those have to always use 1 to code the CENSoring event).
#' Note that factors are not allowed in \code{CENS}.
#' @param TRT A column name in \code{data} for the exposure/treatment variable(s).
#' @param MONITOR A column name in \code{data} for the indicator(s) of monitoring events.
#' @param OUTCOME A column name in \code{data} for the survival OUTCOME variable name, code as 1 for the outcome event.
#' @param noCENScat The level (integer) that indicates CONTINUATION OF FOLLOW-UP for ALL censoring variables. Defaults is 0.
#' Use this to modify the default reference category (no CENSoring / continuation of follow-up)
#' for variables specifed in \code{CENS}.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#' @return ...
# @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
importData <- function(data, ID = "Subject_ID", t_name = "time_period", covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y", noCENScat = 0L, verbose = getOption("stremr.verbose")) {
  gvars$verbose <- verbose
  gvars$noCENScat <- noCENScat
  if (verbose) {
    current.options <- capture.output(str(gvars$opts))
    print("Using the following stremr options/settings: ")
    cat('\n')
    cat(paste0(current.options, collapse = '\n'), '\n')
  }
  if (missing(covars)) { # define time-varing covars (L) as everything else in data besides these vars
    covars <- setdiff(colnames(data), c(ID, t_name, CENS, TRT, MONITOR, OUTCOME))
  }
  # The ordering of variables in this list is the assumed temporal order!
  nodes <- list(Lnodes = covars, Cnodes = CENS, Anodes = TRT, Nnodes = MONITOR, Ynode = OUTCOME, IDnode = ID, tnode = t_name)
  OData <- DataStorageClass$new(Odata = data, nodes = nodes, noCENScat = noCENScat)

  # --------------------------------------------------------------------------------------------------------
  # Check no extra rows after event:
  # --------------------------------------------------------------------------------------------------------
  OData$check_norows_after_event()


  # --------------------------------------------------------------------------------------------------------
  # Create dummies for factors
  # --------------------------------------------------------------------------------------------------------
  factor.Ls <- as.character(CheckExistFactors(OData$dat.sVar))
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
  logical.Ls <- unlist(lapply(OData$dat.sVar, is.logical))
  logical.Ls <- names(logical.Ls)[logical.Ls]
  if (length(logical.Ls)>0 && verbose) message("...converting logical columns to binary integers (0 = FALSE)...")
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


# ---------------------------------------------------------------------------------------
#' Define regression models
#'
#' Alternative approach to manually define parts of the propsensity score model with formula/strata and model controls
#' @param OData OData Input data object created by \code{importData} function.
#' @param regforms Regression formula, only main terms are allowed.
#' @param stratify Expression(s) for creating strata(s) (model will be fit separately on each strata)
#' @param params Additional modeling controls (for downstream modeling algorithms)
#' @export
define_single_regression <- function(OData, regforms, stratify = NULL, params = list()) {
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  return(process_regforms(regforms = regforms, stratify.EXPRS = stratify, model_contrl = params,
                          OData = OData, sVar.map = nodes, factor.map = new.factor.names))
}

internal_define_reg <- function(reg_object, regforms, default.reg, stratify.EXPRS, model_contrl, OData, sVar.map, factor.map, censoring) {
  if (!missing(reg_object)) {
    assert_that(is.list(reg_object))
    if (!"ListOfRegressionForms" %in% class(reg_object)) {
      class(reg_object) <- c(class(reg_object), "ListOfRegressionForms")
    }
    if (censoring) {
      for (reg_idx in seq_along(reg_object)) {
       reg_object[[reg_idx]]$censoring <- TRUE
      }
      reg_object <- stratify_by_uncensored(reg_object)
    }
  } else {
    assert_that(is.list(model_contrl))
    reg_object <- process_regforms(regforms = regforms, default.reg = default.reg, stratify.EXPRS = stratify.EXPRS, model_contrl = model_contrl,
                                   OData = OData, sVar.map = sVar.map, factor.map = factor.map, censoring = censoring)
  }
  return(reg_object)
}

# ---------------------------------------------------------------------------------------
#' Define and fit propensity score models.
#'
#' Defines and fits regression models for the propensity scores for censoring, treatment and monitoring.
#'
#' @param OData Input data object created by \code{importData} function.
#' @param gform_CENS ...
#' @param gform_TRT ...
#' @param gform_MONITOR ...
#' @param stratify_CENS ...
#' @param stratify_TRT ...
#' @param stratify_MONITOR ...
#' @param params_CENS ...
#' @param params_TRT ...
#' @param params_MONITOR ...
#' @param reg_CENS ...
#' @param reg_TRT ...
#' @param reg_MONITOR ...
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#' @return ...
# @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
fitPropensity <- function(OData,
                          gform_CENS, gform_TRT, gform_MONITOR,
                          stratify_CENS = NULL, stratify_TRT = NULL, stratify_MONITOR = NULL,
                          params_CENS = list(), params_TRT = list(), params_MONITOR = list(),
                          reg_CENS, reg_TRT, reg_MONITOR,
                          verbose = getOption("stremr.verbose")) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names

  # ------------------------------------------------------------------------------------------------
  # Process the input formulas and stratification settings;
  # Define regression classes for g.C, g.A, g.N and put them in a single list of regressions.
  # ------------------------------------------------------------------------------------------------
  gform_CENS.default <- "Cnodes ~ Lnodes"
  gform_TRT.default <- "Anodes ~ Lnodes"
  gform_MONITOR.default <- "Nnodes ~ Anodes + Lnodes"
  g_CAN_regs_list <- vector(mode = "list", length = 3)
  names(g_CAN_regs_list) <- c("gC", "gA", "gN")
  class(g_CAN_regs_list) <- c(class(g_CAN_regs_list), "ListOfRegressionForms")

  g_CAN_regs_list[["gC"]] <- internal_define_reg(reg_CENS, gform_CENS, default.reg = gform_CENS.default, stratify.EXPRS = stratify_CENS, model_contrl = params_CENS,
                                                 OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = TRUE)
  g_CAN_regs_list[["gA"]] <- internal_define_reg(reg_TRT, gform_TRT, default.reg = gform_TRT.default, stratify.EXPRS = stratify_TRT, model_contrl = params_TRT,
                                                 OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = FALSE)
  g_CAN_regs_list[["gN"]] <- internal_define_reg(reg_MONITOR, gform_MONITOR, default.reg = gform_MONITOR.default, stratify.EXPRS = stratify_MONITOR, model_contrl = params_MONITOR,
                                                 OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = FALSE)

  # ------------------------------------------------------------------------------------------
  # DEFINE a single regression class
  # Perform S3 method dispatch on ALL_g_regs, which will determine the nested tree of SummaryModel objects
  # Perform fit and prediction
  # ------------------------------------------------------------------------------------------
  modelfits.g0 <- GenericModel$new(reg = g_CAN_regs_list, DataStorageClass.g0 = OData)

  # load data into h2o (not used):
  # mainH2Oframe <- OData$fast.load.to.H2O(OData$dat.sVar,
  #                                       saveH2O = TRUE,
  #                                       destination_frame = "H2OMainDataTable")
  modelfits.g0$fit(data = OData, predict = TRUE)

  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0
  h_gN <- try(modelfits.g0$predictAeqa(n = OData$nobs), silent = TRUE)
  if (inherits(h_gN, "try-error")) { # if failed, it means that prediction cannot be done with newdata
    h_gN <- modelfits.g0$predictAeqa(newdata = OData, n = OData$nobs)
  }

  # ------------------------------------------------------------------------------------------
  # Observed likelihood of (A,C,N) at each t, based on fitted object models in object modelfits.g0
  # ------------------------------------------------------------------------------------------
  # get back g_CAN_regs_list:
  OData$modelfits.g0 <- modelfits.g0
  ALL_g_regs <- modelfits.g0$reg

  OData$modelfit.gC <- modelfits.g0$getPsAsW.models()[[which(names(ALL_g_regs) %in% "gC")]]
  OData$modelfit.gA <- modelfits.g0$getPsAsW.models()[[which(names(ALL_g_regs) %in% "gA")]]
  OData$modelfit.gN <- modelfits.g0$getPsAsW.models()[[which(names(ALL_g_regs) %in% "gN")]]

  g0.A <- OData$modelfit.gA$getcumprodAeqa()
  g0.C <- OData$modelfit.gC$getcumprodAeqa()
  g0.N <- OData$modelfit.gN$getcumprodAeqa()

  OData$dat.sVar[, c("g0.A", "g0.C", "g0.N", "g0.CAN") := list(g0.A, g0.C, g0.N, g0.A*g0.C*g0.N)]
  # newdat <- OData$dat.sVar[, list("g0.A" = g0.A, "g0.C" = g0.C, "g0.N" = g0.N, "g0.CAN" = g0.A*g0.C*g0.N)]
  return(OData)
}

defineNodeGstarIPW <- function(OData, intervened_NODE, NodeNames, useonly_t_NODE, g.obs) {
  # probability of P(A^*(t)=n(t)) or P(N^*(t)=n(t)) under counterfactual A^*(t) or N^*(t) and observed a(t) or n(t)
  # if intervened_NODE returns more than one rule-column, evaluate g^* for each and the multiply to get a single joint (for each time point)
  if (!is.null(intervened_NODE)) {
    for (intervened_NODE_col in intervened_NODE) CheckVarNameExists(OData$dat.sVar, intervened_NODE_col)
    assert_that(length(intervened_NODE) == length(NodeNames))
    # From intervened_NODE we need to evaluate the likelihood: g^*(A^*(t)=A(t)) based on the observed data A(t) and counterfactuals A^*(t) in intervened_NODE
    Q_regs_list <- vector(mode = "list", length = length(NodeNames))
    names(Q_regs_list) <- c(NodeNames)
    class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")
    for (i in seq_along(NodeNames)) {
      reg <- RegressionClass$new(outvar = NodeNames[i], predvars = NULL, outvar.class = list("deterministic"),
                                 subset_vars = list(NodeNames[i]),
                                 model_contrl = list(gstar.Name = intervened_NODE[i]))
      Q_regs_list[[i]] <- reg
    }
    gstar.NODE.obj <- GenericModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)
    gstar.NODE <- gstar.NODE.obj$fit(data = OData)$predictAeqa(n = OData$nobs)
    subset_idx <- OData$evalsubst(subset_exprs = useonly_t_NODE)

    if (any(is.na(subset_idx))) {
      stop("the subset index evaluation for the expression '" %+% useonly_t_NODE %+% "' resulted in NAs")
    }

    gstar.NODE[!subset_idx] <- g.obs[!subset_idx]
  } else {
    # use the actual observed exposure probability (no intervention on NODE)
    gstar.NODE <- g.obs
  }
  return(gstar.NODE)
}

# ---------------------------------------------------------------------------------------
#' Inverse Probability Weights.
#'
#' Evaluate the inverse probability weights for up to 3 intervention nodes: \code{CENS}, \code{TRT} and \code{MONITOR}.
#' This is based on the inverse of the propensity score fits for the observed likelihood (g0.C, g0.A, g0.N),
#' multiplied by the indicator of not being censored and the probability of each intervention in \code{intervened_TRT} and \code{intervened_MONITOR}.
#' Requires column name(s) that specify the counterfactual node values or the counterfactual probabilities of each node being 1 (for stochastic interventions).
#' The output is person-specific data with evaluated weights, \code{wts.DT}, only observation-times with non-zero weight are kept
#' Can be one regimen per single run of this block, which are then combined into a list of output datasets with lapply.
#' Alternative is to allow input with several rules/regimens, which are automatically combined into a list of output datasets.
#' @param OData Input data object created by \code{importData} function.
#' @param intervened_TRT Column name in the input data with the probabilities (or indicators) of counterfactual treatment nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on treatment node(s).
#' @param intervened_MONITOR Column name in the input data with probabilities (or indicators) of counterfactual monitoring nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on the monitoring node(s).
#' @param useonly_t_TRT Use for intervening only on some subset of observation and time-specific treatment nodes.
#' Should be a character string with a logical expression that defines the subset of intervention observations.
#' For example, using \code{"TRT == 0"} will intervene only at observations with the value of \code{TRT} being equal to zero.
#' The expression can contain any variable name that was defined in the input dataset.
#' Leave as \code{NULL} when intervening on all observations/time-points.
#' @param useonly_t_MONITOR Same as \code{useonly_t_TRT}, but for monitoring nodes.
#' @param rule_name Optional name for the treatment/monitoring regimen.
#' @return ...
# @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
getIPWeights <- function(OData, intervened_TRT = NULL, intervened_MONITOR = NULL, useonly_t_TRT = NULL, useonly_t_MONITOR = NULL,
                         rule_name = paste0(c(intervened_TRT, intervened_MONITOR), collapse = "")
                         ) {
  getIPWeights_fun_call <- match.call()
  nodes <- OData$nodes
  if (!is.null(useonly_t_TRT)) assert_that(is.character(useonly_t_TRT))
  if (!is.null(useonly_t_MONITOR)) assert_that(is.character(useonly_t_MONITOR))
  # OData$dat.sVar[, c("g0.CAN.compare") := list(h_gN)] # should be identical to g0.CAN

  # ------------------------------------------------------------------------------------------
  # Probabilities of counterfactual interventions under observed (A,C,N) at each t
  # Combine the propensity score for observed (g0.C, g0.A, g0.N) with the propensity scores for interventions (gstar.C, gstar.A, gstar.N):
  # ------------------------------------------------------------------------------------------------------------------------------
  # (1) gstar.CENS: the indicator of not being censored.
  # (2) gstar.TRT: prob of following one treatment rule; and
  # (3) gstar.MONITOR prob following the monitoring regime; and
  # ------------------------------------------------------------------------------------------------------------------------------
  if (is.null(OData$modelfits.g0)) stop("...cannot locate propensity scores in 'OData' object - must run fitPropensity(...) prior to calling this function")
  if (any(!(c("g0.A", "g0.C", "g0.N", "g0.CAN") %in% names(OData$dat.sVar)))) stop("... fatal error; propensity scores were not found in the input dataset, please re-run fitPropensity(...)")


  # indicator that the person is uncensored at each t (continuation of follow-up)
  gstar.CENS = as.integer(OData$eval_uncensored())
  # Likelihood P(A^*(t)=A(t)) under counterfactual intervention A^*(t) on A(t)
  gstar.TRT <- defineNodeGstarIPW(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, OData$dat.sVar[["g0.A"]])
  # Likelihood for monitoring P(N^*(t)=N(t)) under counterfactual intervention N^*(t) on N(t):
  gstar.MONITOR <- defineNodeGstarIPW(OData, intervened_MONITOR, nodes$Nnodes, useonly_t_MONITOR, OData$dat.sVar[["g0.N"]])
  # Save all likelihoods relating to propensity scores in separate dataset:
  wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, nodes$Ynode, "g0.A", "g0.C", "g0.N", "g0.CAN"), with = FALSE] # [wt.by.t > 0, ]
  setkeyv(wts.DT, cols = c(nodes$IDnode, nodes$tnode))
  wts.DT[, "gstar.C" := gstar.CENS]
  wts.DT[, "gstar.A" := gstar.TRT]
  wts.DT[, "gstar.N" := gstar.MONITOR]
  # Joint likelihoood for all 3 node types:
  wts.DT[, "gstar.CAN" := gstar.CENS * gstar.TRT * gstar.MONITOR]
  # Weights by time and cumulative weights by time:
  wts.DT[, "wt.by.t" := gstar.CAN / g0.CAN, by = eval(nodes$IDnode)][, "cum.IPAW" := cumprod(wt.by.t), by = eval(nodes$IDnode)]

  # -------------------------------------------------------------------------------------------
  # Weight stabilization - get emp P(followed rule at time t | followed rule up to now)
  # -------------------------------------------------------------------------------------------
  nIDs <- OData$nuniqueIDs
  # THE ENUMERATOR: the total sum of subjects followed the rule gstar.A at t
  # THE DENOMINATOR: divide above by the total number of subjects who were still at risk of NOT FOLLOWING the rule at t
  # i.e., followed rule at t-1, assume at the first time-point EVERYONE was following the rule (so denominator = n)
  # (The total sum of all subjects who WERE AT RISK at t)
  # (FASTER) Version outside data.table, then merge back results:
  wts.DT[, "rule.follower.gCAN" := as.integer(cumprod(gstar.CAN) > 0), by = eval(nodes$IDnode)]
  n.follow.rule.t <- wts.DT[, list(N.follow.rule = sum(rule.follower.gCAN, na.rm = TRUE)), by = eval(nodes$tnode)]
  wts.DT[, "rule.follower.gCAN" := NULL]
  n.follow.rule.t[, N.risk := shift(N.follow.rule, fill = nIDs, type = "lag")][, stab.P := ifelse(N.risk > 0, N.follow.rule / N.risk, 1)][, cum.stab.P := cumprod(stab.P)]
  n.follow.rule.t[, c("N.risk", "stab.P") := list(NULL, NULL)]
  setkeyv(n.follow.rule.t, cols = nodes$tnode)
  wts.DT <- wts.DT[n.follow.rule.t, on = nodes$tnode]
  setkeyv(wts.DT, cols = c(nodes$IDnode, nodes$tnode))

  # Multiply the weight by stabilization factor (numerator) (doesn't change the estimand for saturated MSMs, since cum.stab.P cancels out, but makes weights smaller):
  # if (stabilize) wts.DT[, "cum.IPAW" := cum.stab.P * cum.IPAW]
  # # Add the observation-specific weights to the weighted outcome, merge in by ID & t
  # if (!is.null(weights)) {
  #   if (!is.data.table(weights) || (length(names(weights))>3) || !all(c(nodes$IDnode,nodes$tnode) %in% names(weights))) {
  #     stop("input 'weights' must be a data.table with 3 columns, two of which must be named as: '" %+%  nodes$IDnode %+% "' and '" %+% nodes$tnode %+% "'.")
  #   }
  #   wt_col_name <- names(weights)[which(!(names(weights) %in% c(nodes$IDnode,nodes$tnode)))[1]]
  #   weights_used <- weights
  #   setkeyv(weights_used, cols = c(nodes$IDnode, nodes$tnode))
  #   wts.DT <- merge(wts.DT, weights_used, all.x = TRUE)
  # }
  # # Multiply the outcome by the current (cumulative) weight cum.IPAW:
  # wts.DT[, "Wt.OUTCOME" := get(nodes$Ynode)*cum.IPAW]
  # # Multiply the outcome by additional user-supplied weights:
  # if (!is.null(weights)) wts.DT[, "Wt.OUTCOME" := Wt.OUTCOME*get(wt_col_name)]

  wts.DT[, "rule.name" := eval(as.character(rule_name))]

  attributes(wts.DT)[['getIPWeights_fun_call']] <- getIPWeights_fun_call
  attributes(wts.DT)[['intervened_TRT']] <- intervened_TRT
  attributes(wts.DT)[['intervened_MONITOR']] <- intervened_MONITOR
  # attributes(wts.DT)[['stabilize']] <- stabilize
  return(wts.DT)
}

process_opt_wts <- function(wts_data, weights, nodes, adjust_outcome = TRUE) {
  if (!is.null(weights)) {
    if (!is.data.table(weights) || (length(names(weights)) > 3) || !all(c(nodes$IDnode,nodes$tnode) %in% names(weights))) {
      stop("input 'weights' must be a data.table with 3 columns, two of which must be named as: '" %+%  nodes$IDnode %+% "' and '" %+% nodes$tnode %+% "'.")
    }
    wt_col_name <- names(weights)[which(!(names(weights) %in% c(nodes$IDnode,nodes$tnode)))[1]]
    setkeyv(weights, cols = c(nodes$IDnode, nodes$tnode))
    wts_data <- merge(wts_data, weights, all.x = TRUE)
    # Multiply the outcome by additional user-supplied weights:
    if ("Wt.OUTCOME" %in% names(wts_data) && adjust_outcome) {
      wts_data[, "Wt.OUTCOME" := Wt.OUTCOME * get(wt_col_name)]
    } else if (!adjust_outcome) {
      wts_data[, "cum.IPAW" := cum.IPAW * get(wt_col_name)]
    }
  }
  return(wts_data)
}

# ---------------------------------------------------------------------------------------
#' Direct (bounded) IPW estimator of discrete survival function.
#' @param wts_data \code{data.table} returned by a single call to \code{getIPWeights}. Must contain the treatment/monitoring estimated IPTW weights for a SINGLE rule.
#' @param OData The object returned by function \code{fitPropensity}. Contains the input data and the previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param weights (NOT IMPLEMENTED) Optional \code{data.table} with additional observation-time-specific weights.  Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param trunc_weights (NOT IMPLEMENTED) Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @return A data.table with bounded IPW estimates of risk and survival by time.
#' @example tests/examples/2_building_blocks_example.R
#' @export
survDirectIPW <- function(wts_data, OData, weights, trunc_weights) {
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  ## Extract relevant information
  ID.t.IPW.Y <- wts_data[,list(get(nodes$IDnode), get(t_name), cum.IPAW, get(Ynode))]
  names(ID.t.IPW.Y) <- c(nodes$IDnode, t_name, "cum.IPAW", Ynode)

  ## Make sure every patient has an entry for every time point by LVCF
  tmax <- wts_data[, max(get(t_name))]

  tmax <- tmax - 1
  UID <- wts_data[, unique(get(nodes$IDnode))]

  all.ID.t <- as.data.table(cbind(rep(UID,each=(tmax+1)), rep(0:tmax,length(UID)) ))
  names(all.ID.t) <- c(nodes$IDnode, t_name)

  all.ID.t <- merge(all.ID.t, ID.t.IPW.Y, all.x=TRUE, by = c(nodes$IDnode, t_name))
  all.ID.t[ , c("cum.IPAW", Ynode) := list(zoo::na.locf(cum.IPAW), zoo::na.locf(get(Ynode))), by = get(nodes$IDnode)]

  ## Numerator of bounded IPW for survival:
  numIPW <- all.ID.t[, {sum_Y_IPAW = sum(get(Ynode)*cum.IPAW, na.rm = TRUE); list(sum_Y_IPAW = sum_Y_IPAW)}, by = eval(t_name)]
  # numIPW <- all.ID.t[, .(sum_Y_IPAW = sum(get(Ynode)*eval(as.name("cum.IPAW")), na.rm = TRUE)), by = eval(t_name)]

  ## Denominator of bounded IPW for survival:
  denomIPW <- all.ID.t[, {sum_IPAW = sum(cum.IPAW, na.rm = TRUE); list(sum_IPAW = sum_IPAW)}, by = eval(t_name)]
  # denomIPW <- all.ID.t[, .(sum_IPAW = sum(eval(as.name("cum.IPAW")), na.rm = TRUE)), by = eval(t_name)]

  ## Bounded IPW of survival (direct):
  risk.t <- (numIPW[, "sum_Y_IPAW", with = FALSE] / denomIPW[, "sum_IPAW", with = FALSE])
  # S.t.n <- 1 - (numIPW[, "sum_Y_IPAW", with = FALSE] / denomIPW[, "sum_IPAW", with = FALSE])

  resultDT <- data.table(est_name = "DirectBoundedIPW", merge(numIPW, denomIPW, by = t_name))
  resultDT[, c("risk", "S.t.n") := list(risk.t[[1]], 1 - risk.t[[1]])]
  return(resultDT)
}

# ---------------------------------------------------------------------------------------
#' Non-parametric (saturated) MSM for survival based on previously evaluated IP weights.
#' @param wts_data \code{data.table} returned by a single call to \code{getIPWeights}. Must contain the treatment/monitoring estimated IPTW weights for a SINGLE rule.
#' @param OData The object returned by function \code{fitPropensity}. Contains the input data and the previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param weights Optional \code{data.table} with additional observation-time-specific weights.  Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @return A data.table with hazard and survival function estimates by time. Also include the unadjusted Kaplan-Maier estimates.
#' @example tests/examples/2_building_blocks_example.R
#' @export
survNPMSM <- function(wts_data, OData, weights = NULL, trunc_weights = 10^6) {
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode
  rule.name <- unique(wts_data[["rule.name"]])
  if (length(rule.name)>1) stop("wts_data must contain the weights for a single rule, found more than one unique rule name under in 'rule.name' column")

  wts_data_used <- wts_data[,c(nodes$IDnode,nodes$tnode,nodes$Ynode,"cum.stab.P","cum.IPAW"), with = FALSE]
  setkeyv(wts_data_used, cols = c(nodes$IDnode, nodes$tnode))

  # Initialize weighted outcome 'Wt.OUTCOME' to Ynode:
  wts_data_used[, "Wt.OUTCOME" := get(nodes$Ynode)]

  # Add the observation-specific weights to the weighted outcome, merge in by ID & t:
  wts_data_used <- process_opt_wts(wts_data_used, weights, nodes)

  # ------------------------------------------------------------------------
  # CRUDE HAZARD AND KM ESTIMATE OF SURVIVAL:
  # ------------------------------------------------------------------------
  ht.crude <- wts_data_used[cum.IPAW > 0, {ht.KM = sum(Wt.OUTCOME, na.rm = TRUE) / .N; list(ht.KM = ht.KM)}, by = eval(t_name)][, St.KM := cumprod(1 - ht.KM)]
  # ht.crude2 <- wts_data_used[cum.IPAW > 0, .(ht.KM = sum(eval(as.name("Wt.OUTCOME")), na.rm = TRUE) / .N), by = eval(t_name)][, St.KM := cumprod(1 - ht.KM)]
  # ht.crude <- wts_data_used[cum.IPAW > 0, .(ht.KM = sum(eval(as.name(Ynode)), na.rm = TRUE) / .N), by = eval(t_name)][, St.KM := cumprod(1 - ht.KM)]
  setkeyv(ht.crude, cols = t_name)

  # ------------------------------------------------------------------------
  # IPW-ADJUSTED KM (SATURATED MSM):
  # ------------------------------------------------------------------------
  # Multiply the weight by stabilization factor (numerator) (doesn't change the estimand in saturated MSMs, but makes weights smaller):
  wts_data_used[, "cum.IPAW" := cum.stab.P * cum.IPAW]

  # If trunc_weights < Inf then truncate the weights:
  if (trunc_weights < Inf) wts_data_used[cum.IPAW > trunc_weights, cum.IPAW := trunc_weights]

  # Multiply the outcome by cumulative weights in cum.IPAW:
  wts_data_used[, "Wt.OUTCOME" := Wt.OUTCOME * cum.IPAW]
  # wts_data_used[, "Wt.OUTCOME" := get(nodes$Ynode)*cum.IPAW]

  # THE ENUMERATOR FOR THE HAZARD AT t: the weighted sum of subjects who had experienced the event at t:
  sum_Ywt <- wts_data_used[, {sum_Y_IPAW = sum(Wt.OUTCOME, na.rm = TRUE); list(sum_Y_IPAW = sum_Y_IPAW)}, by = eval(t_name)]; setkeyv(sum_Ywt, cols = t_name)
  # sum_Ywt <- wts_data_used[, .(sum_Y_IPAW = sum(Wt.OUTCOME, na.rm = TRUE)), by = eval(t_name)]; setkeyv(sum_Ywt, cols = t_name)

  # THE DENOMINATOR FOR THE HAZARD AT t: The weighted sum of all subjects who WERE AT RISK at t (equivalent to summing cumulative weights cum.IPAW by t):
  sum_Allwt <- wts_data_used[, {sum_all_IPAW = sum(cum.IPAW, na.rm = TRUE); list(sum_all_IPAW = sum_all_IPAW)}, by = eval(t_name)]; setkeyv(sum_Allwt, cols = t_name)
  # sum_Allwt <- wts_data_used[, .(sum_all_IPAW = sum(cum.IPAW, na.rm = TRUE)), by = eval(t_name)]; setkeyv(sum_Allwt, cols = t_name)

  # EVALUATE THE DISCRETE HAZARD ht AND SURVIVAL St OVER t
  St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, ("St.IPTW") := cumprod(1 - ht)]
  # St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("St.IPTW") := .(cumprod(1 - ht))]

  St_ht_IPAW <- merge(St_ht_IPAW, ht.crude, all=TRUE)
  St_ht_IPAW[, "rule.name" := rule.name]
  return(list(wts_data = wts_data_used, trunc_weights = trunc_weights, IPW_estimates = data.frame(St_ht_IPAW)))
}

format_wts_data <- function(wts_data) {
  # Stack the weighted data sets, if those came in a list:
  if (is.data.table(wts_data)) {
    # ...do nothing...
  } else if (is.list(wts_data)) {
    assert_that(all(sapply(wts_data, is.data.table)))
    wts_data <- rbindlist(wts_data)
  } else {
    stop("...wts_data arg must be either a list of rule-specific weight data.tables or one combined weight data.table...")
  }
return(wts_data)
}

# ---------------------------------------------------------------------------------------
#' Estimate Survival with a particular MSM for the survival-hazard function using previously fitted weights.
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the treatment/exposure and monitoring processes based on
#' the user-specified Marginal Structural Model (MSM) for the counterfactual survival function.
#'
#' This routine will run the weighted logistic regression using the (possibly-weighted) outcomes from many regimens, with dummy indicators for each treatment/monitoring
#' regimen available in \code{wts_data} and each follow-up time interval specified in \code{t_breaks}.
#' When \code{use_weights = TRUE}, the logistic regression for the survival hazard is weighted by the \strong{IPW} (Inverse Probability-Weighted or Horvitz-Thompson) estimated weights
#' in \code{wts_data}. These IPW weights are based on the previously fitted propensity scores (function \code{fitPropensity}), allowing
#' adjustment for confounding by possibly non-random assignment to exposure and monitoring and possibly informative right-censoring.
#' @param wts_data A list of \code{data.table}s, each data set is a result of calling the function \code{getIPWeights}. Must contain the treatment/monitoring rule-specific estimated IPTW weights.
#' This argument can be also a single \code{data.table} obtained with \code{data.table::rbindlist(wts_data)}.
#' @param OData The object returned by function \code{fitPropensity}. Contains the input dat and the previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param t_breaks The vector of integer (or numeric) breaks that defines the dummy indicators of the follow-up in the observed data. Used for fitting the parametric (or saturated) MSM for the survival hazard. See Details.
#' @param use_weights Logical value. Set to \code{FALSE} to ignore the weights in \code{wts_data} and fit a "crude" MSM that does not adjust for the possible confounding due to non-random
#' assignment of the exposure/censoring and monitoring indicators.
#' @param stabilize Set to \code{TRUE} for weight stabilization
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @param weights Optional \code{data.table} with additional observation-time-specific weights.  Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param getSEs A logical indicator. Set to \code{TRUE} to evaluate the standard errors for the estimated survival by using the MSM influence curve.
#' @param est_name A string naming the current MSM estimator. Ignored by the current routine but is used when generating reports with \code{make_report_rmd}.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#'
#' @section Details:
#' **********************************************************************
#'
#' \code{t_breaks} is used for defining the time-intervals of the MSM coefficients for estimation of the survival hazard function.
#' The first value in \code{t_breaks} defines a dummy variable (indicator) for a fully closed interval, with each subsequent value in \code{t_breaks} defining a single right-closed time-interval.
#' For example, \code{t_breaks = c(0,1)} will define the MSM dummy indicators: I(min(t) <= t <=0 ), I(0 < t <= 1) and I(1 < t <= max(t)).
#' On the other hand \code{t_breaks = c(1)} will define the following (more parametric) MSM dummy indicators: I(min(t) <= t <=1 ) and I(1 < t <= max(t)).
#' If omitted, the default is to define a saturated (non-parametric) MSM with a separate dummy variable for every unique period in the observed data.
#'
#' @section MSM for the hazard:
#' **********************************************************************
#'
#' @return A named list with items containing the MSM estimation results:
#'  \itemize{
#'  \item \code{est_name} - .
#'  \item \code{St} - .
#'  \item \code{ht} - .
#'  \item \code{MSM.fit} - .
#'  \item \code{MSM.intervals} - .
#'  \item \code{IC.Var.S.d} - .
#'  \item \code{nID} - .
#'  \item \code{wts_data} - .
#'  \item \code{use_weights} - .
#'  \item \code{trunc_weights} - .
#' }
#' @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/4_survMSM_example.R
#' @export
survMSM <- function(wts_data, OData, t_breaks, use_weights = TRUE, stabilize = TRUE, trunc_weights = 10^6, weights = NULL, getSEs = TRUE, est_name = "IPW", verbose = getOption("stremr.verbose")) {
  gvars$verbose <- verbose
  nID <- OData$nuniqueIDs
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  wts_data <- format_wts_data(wts_data)
  rules_TRT <- sort(unique(wts_data[["rule.name"]]))

  if (verbose) print("performing estimation for the following TRT/MONITOR rules found in column 'rule.name': " %+% paste(rules_TRT, collapse=","))

  # Remove all observations with 0 cumulative weights & copy the weights data.table
  wts_data_used <- wts_data[!is.na(cum.IPAW) & !is.na(eval(as.name(Ynode))) & (cum.IPAW > 0), ]
  setkeyv(wts_data_used, cols = c(nodes$IDnode, nodes$tnode))

  # Multiply the weight by stabilization factor (numerator) (doesn't do anything for saturated MSMs, since cum.stab.P cancels out):
  if (stabilize) wts_data_used[, "cum.IPAW" := cum.stab.P * cum.IPAW]

  # If trunc_weights < Inf, do truncation of the weights
  if (trunc_weights < Inf) wts_data_used[cum.IPAW > trunc_weights, cum.IPAW := trunc_weights]

  # Add additional (user-supplied) observation-specific weights to the cumulative weights:
  wts_data_used <- process_opt_wts(wts_data_used, weights, nodes, adjust_outcome = FALSE)

  # When !use_weights run a crude estimator by setting all non-zero weights to 1
  if (!use_weights) wts_data_used[cum.IPAW > 0, cum.IPAW := 1L]

  # Define all observed sequence of periods (t's)
  mint <- min(wts_data_used[[t_name]], na.rm = TRUE); maxt <- max(wts_data_used[[t_name]], na.rm = TRUE)
  periods <- (mint:maxt)
  periods_idx <- seq_along(periods)
  if (verbose) {
    print("periods"); print(periods)
  }

  # Default t_breaks, error checks for t_breaks, plus padding w/ mint & maxt:
  if (missing(t_breaks)) {
    # default t_breaks is to use a saturated (non-parametric) MSM
    t_breaks <- sort(periods)
    if (verbose)
      message("running with default 't_breaks': (" %+%
        paste0(t_breaks, collapse = ",") %+%
        "); \nNote: such 't_breaks' define a separate coefficient for every unique follow-up time period resulting in a saturated (non-parametric) MSM.")
  }
  if (length(unique(t_breaks)) < length(t_breaks)) {
    stop("all t_breaks must be unique")
  }
  if (!all(t_breaks %in% periods)) {
    stop("all t_breaks must be contained between minimum and maximum follow-up periods:" %+% t_breaks[!(t_breaks %in% periods)])
  }

  if (max(t_breaks) < maxt) t_breaks <- sort(c(t_breaks, maxt)) # pad on the right (if needed with maxt):

  # Create the dummies I(d == intervened_TRT) for the logistic MSM for d-specific hazard
  all.d.dummies <- NULL
  for( dummy.j in rules_TRT ){
    wts_data_used[, (dummy.j) := as.integer(rule.name %in% dummy.j)]
    all.d.dummies <- c(all.d.dummies, dummy.j)
  }

  # Create the dummies I(t in interval.j), where interval.j defined by intervals of time of increasing length
  all.t.dummies <- NULL
  t_breaks.mint <- c(mint, t_breaks) # pad t_breaks on the left (with mint)
  MSM.intervals <- matrix(NA, ncol = 2, nrow = length(t_breaks)) # save the actual intervals
  colnames(MSM.intervals) <- c("min.t", "max.t")
  t.per.inteval <- vector(mode = "list", length = nrow(MSM.intervals)) # save the vector of period vals that belong to each interval
  for (t.idx in 2:length(t_breaks.mint)) {
    low.t <- t_breaks.mint[t.idx - 1]
    high.t <- t_breaks.mint[t.idx]
    # First interval needs to be closed on both sides (includes the smallest obesrved follow-up, mint)
    if (t.idx == 2L) {
      dummy.j <- paste("Periods.", low.t, "to", high.t, sep="")
      MSM.intervals[t.idx - 1, ] <- c(low.t, high.t); t.per.inteval[[t.idx - 1]] <- unique(low.t:high.t)
      wts_data_used[, (dummy.j) := as.integer(get(t_name) >= low.t & get(t_name) <= high.t)]
      # wts_data_used[, (dummy.j) := as.integer(eval(as.name(t_name)) >= low.t & eval(as.name(t_name)) <= high.t)]
    } else {
      dummy.j <- paste("Periods.", (low.t + 1), "to", high.t, sep="")
      MSM.intervals[t.idx - 1, ] <- c(low.t + 1, high.t)
      t.per.inteval[[t.idx - 1]] <- unique((low.t+1):high.t)
      wts_data_used[, (dummy.j) := as.integer(get(t_name) >= (low.t + 1) & get(t_name) <= high.t)]
      # wts_data_used[, (dummy.j) := as.integer(eval(as.name(t_name)) >= (low.t + 1) & eval(as.name(t_name)) <= high.t)]
    }
    print("defined t.dummy: " %+% dummy.j)
    all.t.dummies <- c(all.t.dummies, dummy.j)
  }

  # Create interaction dummies I(t in interval.j & d == intervened_TRT)
  for (d.dummy in all.d.dummies) {
    for (t.dummy in all.t.dummies) {
      if (verbose)
        print(t.dummy %+% "_" %+% d.dummy)
      wts_data_used[, (t.dummy %+% "_" %+% d.dummy) := as.integer(eval(as.name(t.dummy)) & eval(as.name(d.dummy)))]
    }
  }

  all_dummies <-  paste(sapply(all.d.dummies, function(x) {
                        return(paste(paste(paste(all.t.dummies, x, sep="_"), sep="")))
                        }))

  # Fit the hazard MSM
  resglmMSM <- runglmMSM(OData, wts_data_used, all_dummies, Ynode, verbose)
  wts_data_used[, glm.IPAW.predictP1 := resglmMSM$glm.IPAW.predictP1]
  m.fit <- resglmMSM$m.fit

  # Compute the Survival curves under each d
  S2.IPAW <- hazard.IPAW <- rep(list(rep(NA,maxt-mint+1)), length(rules_TRT))
  names(S2.IPAW) <- names(hazard.IPAW) <- rules_TRT

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
      # print("fetching MSM coefficient for period " %+% period.j %+% " and rule " %+% d.j %+% ": " %+% MSM.term)
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
    design.d.t <- rep(list(matrix(0L, ncol = length(all_dummies), nrow = length(mint:maxt))), length(rules_TRT))
    IC.Var.S.d <- vector(mode = "list", length(rules_TRT))
    names(design.d.t) <- names(IC.Var.S.d) <- rules_TRT
    # the matrix where each row consists of indicators for t-specific derivatives of m(t,d), for each fixed d.
    # the rows loop over all possible t's for which the survival will be plotted! Even if there was the same coefficient beta for several t's
    # p.coef - number of time-specific coefficients in the MSM
    p.coef <- nrow(MSM.intervals) # p.coef <- length(tjmin)
    design.t <- matrix(0L, ncol = p.coef, nrow = length(periods))
    for (period.idx in seq_along(periods)) {
      period.j <- periods[period.idx]
      col.idx <- which(period.j <= MSM.intervals[,2] & period.j >= MSM.intervals[,1])
      design.t[period.idx, col.idx] <- 1
    }
    beta.IC.O.SEs <- getSEcoef(ID = nodes$IDnode, nID = nID, t.var = nodes$tnode, Yname = Ynode,
                              MSMdata = wts_data_used, MSMdesign = as.matrix(wts_data_used[, all_dummies, with = FALSE]),
                              MSMpredict = "glm.IPAW.predictP1", IPW_MSMestimator = use_weights)

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
              est_name = est_name,
              periods = periods,
              St = S2.IPAW,
              ht = hazard.IPAW,
              MSM.fit = m.fit,
              MSM.intervals = MSM.intervals,
              IC.Var.S.d = IC.Var.S.d,
              nID = nID,
              nobs = nrow(wts_data_used),
              wts_data = wts_data_used,
              use_weights = use_weights,
              trunc_weights = trunc_weights
            )
  return(MSM_out)
}

runglmMSM <- function(OData, wts_data, all_dummies, Ynode, verbose) {
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
      MSM.designmat.H2O <- OData$fast.load.to.H2O(wts_data,
                                                  saveH2O = FALSE,
                                                  destination_frame = "MSM.designmat.H2O")
    )
    if (verbose) { print("time to load the design mat into H2OFRAME: "); print(loadframe_t) }

    m.fit_h2o <- try(h2o::h2o.glm(y = Ynode,
                                  x = all_dummies,
                                  intercept = FALSE,
                                  weights_column = "cum.IPAW",
                                  training_frame = MSM.designmat.H2O,
                                  family = "binomial",
                                  standardize = FALSE,
                                  solver = c("IRLSM"), # solver = c("L_BFGS"),
                                  lambda = 0L,
                                  max_iterations = 50,
                                  ignore_const_cols = FALSE
                                  ),
                silent = TRUE)

    out_coef <- vector(mode = "numeric", length = length(all_dummies))
    out_coef[] <- NA
    names(out_coef) <- c(all_dummies)
    out_coef[names(m.fit_h2o@model$coefficients)[-1]] <- m.fit_h2o@model$coefficients[-1]
    m.fit <- list(coef = out_coef, linkfun = "logit_linkinv", fitfunname = "h2o.glm")
    glm.IPAW.predictP1 <- as.vector(h2o::h2o.predict(m.fit_h2o, newdata = MSM.designmat.H2O)[,"p1"])
    # wts_data[, glm.IPAW.predictP1 := as.vector(h2o::h2o.predict(m.fit_h2o, newdata = MSM.designmat.H2O)[,"p1"])]
  } else {
    if (verbose) message("...fitting hazard MSM with speedglm::speedglm.wfit...")
    Xdesign.mat <- as.matrix(wts_data[, all_dummies, with = FALSE])
    m.fit <- try(speedglm::speedglm.wfit(
                                       X = Xdesign.mat,
                                       y = as.numeric(wts_data[[Ynode]]),
                                       intercept = FALSE,
                                       family = quasibinomial(),
                                       weights = wts_data[["cum.IPAW"]],
                                       trace = FALSE),
                        silent = TRUE)
    if (inherits(m.fit, "try-error")) { # if failed, fall back on stats::glm
      if (verbose) message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", m.fit)
      ctrl <- glm.control(trace = FALSE)
      SuppressGivenWarnings({
        m.fit <- stats::glm.fit(x = Xdesign.mat,
                                y = as.numeric(wts_data[[Ynode]]),
                                family = quasibinomial(),
                                intercept = FALSE, control = ctrl)
      }, GetWarningsToSuppress())
    }
    m.fit <- list(coef = m.fit$coef, linkfun = "logit_linkinv", fitfunname = "speedglm")
    if (verbose) {
      print("MSM fits"); print(m.fit$coef)
    }
    glm.IPAW.predictP1 <- logispredict(m.fit, Xdesign.mat)
    # wts_data[, glm.IPAW.predictP1 := logispredict(m.fit, Xdesign.mat)]
  }
  return(list(glm.IPAW.predictP1 = glm.IPAW.predictP1, m.fit = m.fit))
  # return(list(wts_data = wts_data, m.fit = m.fit))
}