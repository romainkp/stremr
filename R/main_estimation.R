# adding to appease CRAN check with non-standard eval in data.table:

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("gstar.CAN", "g0.CAN", "wt.by.t", "rule.follower.gCAN", "new.TRT.gstar",
                            "N.risk", "N.follow.rule", "stab.P", "cum.stab.P", "cum.IPAW",
                            "rule.name", "glm.IPAW.predictP1", "St.KM", "Wt.OUTCOME",
                            "St_ht_IPAW", "ht.NPMSM", "ht.KM", "EIC_i_t0", "EIC_i_tplus",
                            "g0.A", "g0.C", "g0.N"))
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
    reg_object <- process_regforms(regforms = regforms,
                                   default.reg = default.reg,
                                   stratify.EXPRS = stratify.EXPRS,
                                   model_contrl = model_contrl,
                                   OData = OData,
                                   sVar.map = sVar.map,
                                   factor.map = factor.map,
                                   censoring = censoring)
  }
  return(reg_object)
}

# ---------------------------------------------------------------------------------------
#' Directly specify a single regression model
#'
#' This function provides an alternative way to manually define parts of the propsensity score model
#' with formula/strata and model controls.
#' This function is for advanced users. It provides explicit and manual control
#' over every single model fit, e.g., every strata of the exposure propensity score.
#' @param OData OData Input data object created by \code{importData} function.
#' @param regforms Regression formula, only main terms are allowed.
#' @param stratify Expression(s) for creating strata(s) (model will be fit separately on each strata)
#' @param models Optional parameter specifying the models with \code{gridisl} R package.
#' Must be an object of class \code{ModelStack} specified with \code{gridisl::defModel} function.
#' @param fit_method Model selection approach. Can be \code{"none"} - no model selection,
#' \code{"cv"} - V fold cross-validation that selects the best model according to lowest cross-validated MSE
#' (must specify the column name that contains the fold IDs).
# \code{"holdout"} - model selection by splitting the data into training and validation samples according to
# lowest validation sample MSE (must specify the column of \code{TRUE} / \code{FALSE} indicators,
# where \code{TRUE} indicates that this row will be selected as part of the model validation sample).
#' @param fold_column The column name in the input data (ordered factor) that contains the fold IDs to be used as part of the validation sample.
#' Use the provided function \code{\link{define_CVfolds}} to
#' define such folds or define the folds using your own method.
# @param params Additional modeling controls (for downstream modeling algorithms)
#' @param ... Additional parameters that will be passed down directly to the modeling function.
#' @export
define_single_regression <- function(OData,
                                     regforms,
                                     stratify = NULL,
                                     models = NULL,
                                     fit_method = stremrOptions("fit_method"),
                                     fold_column = stremrOptions("fold_column"),
                                     ...) {
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names

  opt_params <- capture.exprs(...)

  if (!is.null(models) && suppressWarnings(!is.na(models))) {
    assert_that(is.ModelStack(models) || is(models, "Lrnr_base"))
  } else {
    models <- do.call(sl3::Lrnr_glm_fast$new, opt_params)
  }

  models_control <- c(list(models = models), 
                           opt_params = list(opt_params)
                      )

  models_control[["fit_method"]] <- fit_method[1L]
  models_control[["fold_column"]] <- fold_column

  return(process_regforms(regforms = regforms,
                          stratify.EXPRS = stratify,
                          model_contrl = models_control,
                          OData = OData,
                          sVar.map = nodes,
                          factor.map = new.factor.names))
}

# ---------------------------------------------------------------------------------------
#' Import data, define various nodes, define dummies for factor columns and define OData R6 object
#'
#' @param data Input data in long format. Can be a \code{data.frame} or a \code{data.table} with named columns,
#' containing the time-varying covariates (\code{covars}),
#' the right-censoring event indicator(s) (\code{CENS}), the exposure variable(s) (\code{TRT}),
#' the monitoring process variable(s) (\code{MONITOR})
#' and the survival OUTCOME variable (\code{OUTCOME}).
#' @param ID Unique subject identifier column name in \code{data}.
#' @param t_name The name of the time/period variable in \code{data}.
#' @param covars Vector of names with time varying and baseline covariates in \code{data}.
#' This argument does not need to be specified, by default all variables
#' that are not in \code{ID}, \code{t}, \code{CENS}, \code{TRT}, \code{MONITOR} and \code{OUTCOME}
#' will be considered as covariates.
#' @param CENS Column name of the censoring variable(s) in \code{data}.
#' Leave as missing if not intervening on censoring / no right-censoring in the input data.
#' Each separate variable specified in \code{CENS} can be either binary (0/1 valued integer) or categorical (integer).
#' For binary indicators of CENSoring, the value of 1 indicates the CENSoring or end of follow-up event
#' (this cannot be changed).
#' For categorical CENSoring variables, by default the value of 0 indicates no CENSoring / continuation of
#' follow-up and other values indicate different reasons for CENSoring.
#' Use the argument \code{noCENScat} to change the reference (continuation of follow-up) category from
#' default 0 to any other value.
#' (NOTE: Changing \code{noCENScat} has zero effect on coding of the binary CENSoring variables, those
#' have to always use 1 to code the CENSoring event).
#' Note that factors are not allowed in \code{CENS}.
#' @param TRT A column name in \code{data} for the exposure/treatment variable(s).
#' @param MONITOR A column name in \code{data} for the indicator(s) of monitoring events.
#' Leave as missing if not intervening on monitoring.
#' @param OUTCOME A column name in \code{data} for the survival OUTCOME variable name, code as 1 for the outcome event.
#' @param weights Optional column name in \code{data} that contains subject- and time-specific weights.
#' @param noCENScat The level (integer) that indicates CONTINUATION OF FOLLOW-UP for ALL censoring variables. Defaults is 0.
#' Use this to modify the default reference category (no CENSoring / continuation of follow-up)
#' for variables specifed in \code{CENS}.
#' @param remove_extra_rows Remove extra rows after the event of interest (survival outcome) has occurred (OUTCOME=1).
#' Set this to FALSE for non-survival data (i.e., when the outcome is not time-to-event and new observations may occur after OUTCOME = 1).
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#' @return ...
# @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
importData <- function(data,
                       ID = "Subject_ID",
                       t_name = "time_period",
                       covars,
                       CENS = NULL, # CENS = "C",
                       TRT = "A",
                       MONITOR = NULL, # MONITOR = "N",
                       OUTCOME = "Y",
                       weights = NULL,
                       noCENScat = 0L,
                       remove_extra_rows = TRUE,
                       verbose = getOption("stremr.verbose")) {
  gvars$verbose <- verbose
  gvars$noCENScat <- noCENScat
  if (verbose) {
    current.options <- capture.output(str(gvars$opts))
    print("stremr will use the following options as defaults: ")
    cat('\n')
    cat(paste0(current.options, collapse = '\n'), '\n')
  }

  if (missing(TRT)) stop("treatment column names must be specified w/ arg 'TRT'")
  if (missing(CENS)) CENS <- NULL
  if (missing(MONITOR)) MONITOR <- NULL

  if (missing(covars)) { # define time-varing covars (L) as everything else in data besides these vars
    covars <- setdiff(colnames(data), c(ID, t_name, CENS, TRT, MONITOR, OUTCOME))
  }
  # The ordering of variables in this list is the assumed temporal order!
  nodes <- list(Lnodes = covars,
                Cnodes = CENS,
                Anodes = TRT,
                Nnodes = MONITOR,
                Ynode = OUTCOME,
                IDnode = ID,
                tnode = t_name,
                weights = weights)
  OData <- DataStorageClass$new(Odata = data, nodes = nodes, noCENScat = noCENScat)

  # --------------------------------------------------------------------------------------------------------
  # Check no extra rows after event:
  # --------------------------------------------------------------------------------------------------------
  if (remove_extra_rows) OData$check_norows_after_event()

  # --------------------------------------------------------------------------------------------------------
  # Create dummies for factors
  # --------------------------------------------------------------------------------------------------------
  factor.Ls <- as.character(CheckExistFactors(OData$dat.sVar))
  new.factor.names <- vector(mode="list", length=length(factor.Ls))
  names(new.factor.names) <- factor.Ls
  if (length(factor.Ls)>0 && verbose)
    message("...converting the following factor(s) to binary dummies (and droping the first factor levels): " %+% paste0(factor.Ls, collapse=","))
  for (factor.varnm in factor.Ls) {
    factor.levs <- levels(OData$dat.sVar[[factor.varnm]])
    ## only define new dummies for factors with > 2 levels
    if (length(factor.levs) > 2) {
      factor.levs <- factor.levs[-1] # remove the first level (reference class)
      # use levels to define cat indicators:
      OData$dat.sVar[,(factor.varnm %+% "_" %+% factor.levs) := lapply(factor.levs, function(x) as.integer(levels(get(factor.varnm))[get(factor.varnm)] %in% x))]
      # to remove the origional factor var: # OData$dat.sVar[,(factor.varnm):=NULL]
      new.factor.names[[factor.varnm]] <- factor.varnm %+% "_" %+% factor.levs
    ## Convert existing factor variable to integer
    } else {
      OData$dat.sVar[, (factor.varnm) := as.integer(levels(get(factor.varnm))[get(factor.varnm)] %in% factor.levs[2])]
      new.factor.names[[factor.varnm]] <- factor.varnm
    }
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

  for (Cnode  in nodes$Cnodes) CheckVarNameExists(OData$dat.sVar, Cnode)
  for (Anode  in nodes$Anodes) CheckVarNameExists(OData$dat.sVar, Anode)
  for (Nnode  in nodes$Nnodes) CheckVarNameExists(OData$dat.sVar, Nnode)
  for (Ynode  in nodes$Ynode)  CheckVarNameExists(OData$dat.sVar, Ynode)
  for (Lnode  in nodes$Lnodes) CheckVarNameExists(OData$dat.sVar, Lnode)
  for (IDnode in nodes$IDnode) CheckVarNameExists(OData$dat.sVar, IDnode)
  for (weights in nodes$weights) CheckVarNameExists(OData$dat.sVar, weights)

  return(OData)
}

define_propensity_model <- function(models, opt_params) {
  ## verify the learners are of acceptible type, otherwise define default learner as sl3 glm
  if (!is.null(models) && suppressWarnings(!is.na(models))) {
    assert_that(is.ModelStack(models) || is(models, "Lrnr_base"))
  } else {
    models <- do.call(sl3::Lrnr_glm_fast$new, opt_params)
  }
  
  models_control <- c(list(models     = models),
                           opt_params = list(opt_params))
  return(models_control)
}

# ---------------------------------------------------------------------------------------
#' Define and fit propensity score models.
#'
#' Defines and fits estimators for the propensity scores, separately for censoring, treatment and monitoring events.
#' When there is right-censoring and/or not intervening on monitoring, only the propensity score model for treatment will be estimated.
#'
#' @param OData Input data object created by \code{importData} function.
#' @param gform_CENS Specify the regression formula for the right-censoring mechanism,
#' in the format "CensVar1 + CensVar2 ~ Predictor1 + Predictor2".
#' Leave as missing for data with no right-censoring.
#' @param gform_TRT Specify the regression formula for the treatment mechanism, in the format "TRTVar1 + TRTVar2 ~ Predictor1 + Predictor2".
#' @param gform_MONITOR  Specify the regression formula for the treatment mechanism, in the format "TRTVar1 + TRTVar2 ~ Predictor1 + Predictor2".
#' Leave as missing for data with no monitoring events or when not intervening on monitoring.
#' @param stratify_CENS Define strata(s) for each censoring variable from \code{gform_CENS}.
#' Must be named list containing the logical expressions (the logical expressions must be provided as character strings).
#' When missing (default), all censoring models will be fit by pooling all available observations, across all time-points.
#' When used the censoring models in \code{gform_CENS} will be trained separately on each strata (defined by separate logical expressions).
#' If the \code{gform_CENS} contains more than one censoring variable
#' then this argumement (\code{stratify_CENS}) must provide separate stratas for each censoring variable or be left as missing.
#' For example, when \code{gform_CENS}="CensVar1 + CensVar2 ~ Predictor1 + Predictor2", this argument should be a list of length two,
#' with list items named as "CensVar1" and "CensVar2". The expressions in stratify_CENS[["CensVar1"]] define the training stratas
#' for censoring variable \code{CensVar1}, while the expressions in stratify_CENS[["CensVar1"]] define the training stratas
#' for \code{CensVar2}. See additional examples below.
#' @param stratify_TRT Define strata(s) for treatment model(s).
#' Must be a list of logical expressions (input the expression as character strings).
#' When missing (default), the treatment model(s) are fit by pooling all available (uncensored) observations, across all time-points.
#' The rules are the same as for \code{stratify_CENS}.
#' @param stratify_MONITOR Define strata(s) for monitoring model(s).
#' Must be a list of logical expressions (input the expression as character strings).
#' When missing (default), the monitoring model is fit by pooling all available (uncensored) observations, across all time-points.
#' The rules are the same as for \code{stratify_CENS}.
#' @param models_CENS Optional parameter specifying the models for fitting the censoring mechanism(s) with
#' \code{gridisl} R package.
#' Must be an object of class \code{ModelStack} specified with \code{gridisl::defModel} function.
#' @param models_TRT Optional parameter specifying the models for fitting the treatment (exposure) mechanism(s)
#' with \code{gridisl} R package.
#' Must be an object of class \code{ModelStack} specified with \code{gridisl::defModel} function.
#' @param models_MONITOR Optional parameter specifying the models for fitting the monitoring mechanism with
#' \code{gridisl} R package.
#' Must be an object of class \code{ModelStack} specified with \code{gridisl::defModel} function.
#' @param fit_method Model selection approach. Can be \code{"none"} - no model selection,
#' \code{"cv"} - V fold cross-validation that selects the best model according to lowest cross-validated MSE
#' (must specify the column name that contains the fold IDs).
# \code{"holdout"} - model selection by splitting the data into training and validation samples according to
# lowest validation sample MSE (must specify the column of \code{TRUE} / \code{FALSE} indicators,
# where \code{TRUE} indicates that this row will be selected as part of the model validation sample).
#' @param fold_column The column name in the input data (ordered factor) that contains the fold IDs to be used as part of the validation sample.
#' Use the provided function \code{\link{define_CVfolds}} to
#' define such folds or define the folds using your own method.
#' @param reg_CENS (ADVANCED FEATURE). Manually define and input the regression specification for each strata of censoring model,
#' using the function \code{\link{define_single_regression}}.
#' @param reg_TRT (ADVANCED FEATURE). Manually define and input the regression specification for each strata of treatment model,
#' using the function \code{\link{define_single_regression}}.
#' @param reg_MONITOR (ADVANCED FEATURE). Manually define and input the regression specification for each strata of monitoring model,
#' using the function \code{\link{define_single_regression}}.
#' @param use_weights (NOT IMPLEMENTED) Set to \code{TRUE} to pass the previously specified weights column to all 
#' learners for all of the propensity score models. This will result in a weights regression models being fit for P(A(t)|L(t)), P(C(t)|L(t)), P(N(t)|L(t)).
#' (NOTE: This will only work when using sl3 learners, such as the default sl3 learner \code{sl3::Lrnr_glm_fast}).
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#' @param ... When all or some of the \code{models_...} arguments are NOT specified, these additional
#' arguments will be passed on directly to all \code{gridisl}
#' modeling functions that are called from this routine,
#' e.g., \code{family = "binomial"} can be used to specify the model family. Note that all such arguments
#' must be named.
#' @return ...
# @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
fitPropensity <- function(OData,
                          gform_CENS,
                          gform_TRT,
                          gform_MONITOR,
                          stratify_CENS = NULL,
                          stratify_TRT = NULL,
                          stratify_MONITOR = NULL,
                          models_CENS = NULL,
                          models_TRT = NULL,
                          models_MONITOR = NULL,
                          fit_method = stremrOptions("fit_method"),
                          fold_column = stremrOptions("fold_column"),
                          reg_CENS,
                          reg_TRT,
                          reg_MONITOR,
                          use_weights = FALSE,
                          # type_CENS = "univariate",
                          # type_TRT = "univariate",
                          # type_MONITOR = "univariate",
                          verbose = getOption("stremr.verbose"),
                          ...) {

  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  opt_params <- capture.exprs(...)
  if (!("family" %in% names(opt_params))) opt_params[["family"]] <- quasibinomial()

  models_CENS_control <-    define_propensity_model(models_CENS,    opt_params)
  models_TRT_control <-     define_propensity_model(models_TRT,     opt_params)
  models_MONITOR_control <- define_propensity_model(models_MONITOR, opt_params)

  models_CENS_control[["fit_method"]] <- models_TRT_control[["fit_method"]] <- models_MONITOR_control[["fit_method"]] <- fit_method[1L]
  models_CENS_control[["fold_column"]] <- models_TRT_control[["fold_column"]] <- models_MONITOR_control[["fold_column"]] <- fold_column

  if (use_weights) stop("...use_weights feature is not implemented yet...")

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

  g_CAN_regs_list[["gC"]] <- internal_define_reg(reg_CENS, gform_CENS,
                                                 default.reg = gform_CENS.default,
                                                 stratify.EXPRS = stratify_CENS,
                                                 model_contrl = models_CENS_control,
                                                 OData = OData,
                                                 sVar.map = nodes,
                                                 factor.map = new.factor.names,
                                                 censoring = TRUE)

  g_CAN_regs_list[["gA"]] <- internal_define_reg(reg_TRT, gform_TRT,
                                                 default.reg = gform_TRT.default,
                                                 stratify.EXPRS = stratify_TRT,
                                                 model_contrl = models_TRT_control,
                                                 OData = OData,
                                                 sVar.map = nodes,
                                                 factor.map = new.factor.names,
                                                 censoring = FALSE)

  g_CAN_regs_list[["gN"]] <- internal_define_reg(reg_MONITOR, gform_MONITOR,
                                                 default.reg = gform_MONITOR.default,
                                                 stratify.EXPRS = stratify_MONITOR,
                                                 model_contrl = models_MONITOR_control,
                                                 OData = OData,
                                                 sVar.map = nodes,
                                                 factor.map = new.factor.names,
                                                 censoring = FALSE)

  # ------------------------------------------------------------------------------------------
  # DEFINE a single regression class
  # Perform S3 method dispatch on ALL_g_regs, which will determine the nested tree of SummaryModel objects
  # Perform fit and prediction
  # ------------------------------------------------------------------------------------------
  modelfits.g0 <- ModelGeneric$new(reg = g_CAN_regs_list, DataStorageClass.g0 = OData)
  modelfits.g0$fit(data = OData, predict = TRUE)

  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0
  h_gN <- try({modelfits.g0$predictAeqa(n = OData$nobs)}, silent = TRUE)

  if (inherits(h_gN, "try-error")) { # if failed, it means that prediction cannot be done without newdata
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

  g_preds <- data.table::data.table(g0.C = OData$modelfit.gC$getcumprodAeqa(),
                                    g0.A = OData$modelfit.gA$getcumprodAeqa(),
                                    g0.N = OData$modelfit.gN$getcumprodAeqa())
  OData$g_preds <- g_preds

  ## warning: side-effect function call, the predicted probabilities are updated inside each class and are saved
  h_gN_holdout <- modelfits.g0$predictAeqa(newdata = OData, n = OData$nobs, holdout = TRUE)
  g_holdout_preds <- data.table::data.table(g0.C = OData$modelfit.gC$getcumprodAeqa(),
                                            g0.A = OData$modelfit.gA$getcumprodAeqa(),
                                            g0.N = OData$modelfit.gN$getcumprodAeqa())
  OData$g_holdout_preds <- g_holdout_preds

  return(OData)
}

## Evaluate the intervention probability  P(A^*(t)=a(t)) / P(N^*(t)=n(t))
## for counterfactual A^*(t) / N^*(t) and the observed data values a(t) / n(t).
## When intervened_NODE contains more than one rule-column, evaluate g^* for each and
## multiply to get a single joint probability (at each time point).
## This function simply grabs the counterfactual node values (N^*(t)) and compares them
## to the observed values (N(t)) by evaluating the indicator (N^*(t)=N(t)).
## The call to fit below is empty, i.e., does nothing other than call ModelDeterministic$predict()
defineNodeGstarIPW <- function(OData, intervened_NODE, NodeNames, useonly_t_NODE, g.obs, modelfit.g, type_intervened) {
  if (!is.null(intervened_NODE) && !is.na(intervened_NODE)) {
    for (intervened_NODE_col in intervened_NODE) CheckVarNameExists(OData$dat.sVar, intervened_NODE_col)
    assert_that(length(intervened_NODE) == length(NodeNames))
    if (!is.null(type_intervened)) {
      assert_that(length(type_intervened)==length(intervened_NODE))
      assert_that(is.character(type_intervened))
      assert_that(all(type_intervened %in% c("bin", "shift", "MSM")))
    } else {
      type_intervened <- rep("bin", length(intervened_NODE))
    }

    # From intervened_NODE we need to evaluate the likelihood: g^*(A^*(t)=A(t)) based on the observed data A(t) and counterfactuals A^*(t) in intervened_NODE
    regs_list <- vector(mode = "list", length = length(NodeNames))
    names(regs_list) <- c(NodeNames)
    class(regs_list) <- c(class(regs_list), "ListOfRegressionForms")
    for (i in seq_along(NodeNames)) {
      modelfit.g.node <- modelfit.g$getPsAsW.models()[[i]]

      reg <- RegressionClass$new(outvar = NodeNames[i],
                                 predvars = NULL,
                                 outvar.class = list("deterministic"),
                                 subset_vars = list(NodeNames[i]),
                                 model_contrl = list(gstar.Name = intervened_NODE[i],
                                                     type_intervened = type_intervened[i],
                                                     modelfit.g = modelfit.g.node))
      regs_list[[i]] <- reg
    }
    gstar.NODE.obj <- ModelGeneric$new(reg = regs_list, DataStorageClass.g0 = OData)
    gstar.NODE <- gstar.NODE.obj$fit(data = OData)$predictAeqa(newdata = OData, n = OData$nobs)
    # gstar.NODE <- gstar.NODE.obj$fit(data = OData)$predictAeqa(n = OData$nobs)
    subset_idx <- OData$evalsubst(subset_exprs = useonly_t_NODE)

    if (any(is.na(subset_idx)))
      stop("the subset index evaluation for the expression '" %+% useonly_t_NODE %+% "' resulted in NAs")

    idx_set_to_g0 <- setdiff(1:length(gstar.NODE), subset_idx)
    gstar.NODE[idx_set_to_g0] <- g.obs[idx_set_to_g0]

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
#' multiplied by the indicator of not being censored and the probability of each intervention in \code{intervened_TRT}
#' and \code{intervened_MONITOR}.
#' Requires column name(s) that specify the counterfactual node values or the counterfactual probabilities of each
#' node being 1 (for stochastic interventions).
#' The output is person-specific data with evaluated weights, \code{wts.DT}, only observation-times with non-zero
#' weight are kept
#' Can be one regimen per single run of this block, which are then combined into a list of output datasets with lapply.
#' Alternative is to allow input with several rules/regimens, which are automatically combined into a list of output datasets.
#' @param OData Input data object created by \code{importData} function.
#' @param intervened_TRT Column name in the input data with the probabilities (or indicators) of counterfactual
#' treatment nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on treatment node(s).
#' @param intervened_MONITOR Column name in the input data with probabilities (or indicators) of counterfactual
#' monitoring nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on the monitoring node(s).
#' @param useonly_t_TRT Use for intervening only on some subset of observation and time-specific treatment nodes.
#' Should be a character string with a logical expression that defines the subset of intervention observations.
#' For example, using \code{"TRT == 0"} will intervene only at observations with the value of \code{TRT} being
#' equal to zero.
#' The expression can contain any variable name that was defined in the input dataset.
#' Leave as \code{NULL} when intervening on all observations/time-points.
#' @param useonly_t_MONITOR Same as \code{useonly_t_TRT}, but for monitoring nodes.
#' @param rule_name Optional name for the treatment/monitoring regimen.
#' @param tmax Maximum value of the follow-up period.
#' All person-time observations above this value will be excluded from the output weights dataset.
#' @param ignore_tmin (ADVANCED FEATURE) Minimum value of the follow-up period at which the IP-weights should start accumulating over time.
#' All IP-weights for time-points t < ignore_tmin will be set to a constant 1.
#' This will have the effect of completely ignoring all weight contributions that occur before ignore_tmin.
#' @param ignore_tmax (ADVANCED FEATURE) Maximum value of the follow-up period to start accumulative the weights over time.
#' All of the time-specific IP-weights with t < ignore_tmin will be set to constant 1 PRIOR to the evaluation of the cumulative weights.
#' This will have the effect of completely ignoring all the IP weight contributions up to and including the time-point ignore_tmin.
#' @param reverse_wt_prod Set to TRUE to take the product of the cumulative weights in reverse time-ordering. That is, the
#' cumulative product will be evaluated by starting from the highest follow-up time point (time variable value).
#' @param holdout Obtain the weights based on out-of-sample (holdout / validation set) predictions of propensity scores.
#' This is useful for running CV-TMLE or evaluating the quality of the model fits based on validation sets.
#' @param eval_stabP Evaluate the additional weight stabilization factor for each time-point.
#' This is used for MSMs only and is enabled by default.
#' @param trunc_weights Specify the numeric weight truncation value. All final weights exceeding the value in
#' \code{trunc_weights} will be truncated.
#' @param type_intervened_TRT (ADVANCED FUNCTIONALITY) Set to \code{NULL} by default, can be characters that are set to either 
#' \code{"bin"}, \code{"shift"} or \code{"MSM"}. 
#' Provides support for different types of interventions on \code{TRT} (treatment) node (counterfactual treatment node \code{A^*(t)}).
#' The default behavior is the same as \code{"bin"}, which assumes that \code{A^*(t)} is binary and 
#' is set equal to either \code{0}, \code{1} or \code{p(t)}, where 0<=\code{p(t)}<=1. 
#' Here, \code{p(t)} denotes the probability that counterfactual A^*(t) is equal to 1, i.e., P(A^*(t)=1)=\code{p(t)} 
#' and it can change in time and subject to subject.
#' For \code{"shift"}, it is assumed that the intervention node \code{A^*(t)} is a shift in the value of the continuous treatment \code{A}, 
#' i.e., \code{A^*(t)}=\code{A(t)}+delta(t).
#' Finally, for "MSM" it is assumed that we simply want the final intervention density \code{g^*(t)} to be set to a constant 1. 
#' This has use for static MSMs. 
#' @param type_intervened_MONITOR (ADVANCED FUNCTIONALITY) Same as \code{type_intervened_TRT}, but for monitoring intervention node 
#' (counterfactual monitoring node \code{N^*(t)}).
#' @return A \code{data.table} with cumulative weights for each subject and each time-point saved under column "cum.IPAW".
# @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/2_building_blocks_example.R
#' @export
getIPWeights <- function(OData,
                         intervened_TRT = NULL,
                         intervened_MONITOR = NULL,
                         useonly_t_TRT = NULL,
                         useonly_t_MONITOR = NULL,
                         rule_name = paste0(c(intervened_TRT, intervened_MONITOR), collapse = ""),
                         tmax = NULL,
                         ignore_tmin = NULL,
                         ignore_tmax = NULL,
                         reverse_wt_prod = FALSE,
                         holdout = FALSE,
                         eval_stabP = TRUE,
                         trunc_weights = Inf,
                         type_intervened_TRT = NULL,
                         type_intervened_MONITOR = NULL
                         ) {
  getIPWeights_fun_call <- match.call()
  nodes <- OData$nodes

  if (!holdout) {
    g_preds <- OData$g_preds
  } else {
    g_preds <- OData$g_holdout_preds
  }

  if (!is.null(useonly_t_TRT) && !is.na(useonly_t_TRT)) assert_that(is.character(useonly_t_TRT))
  if (!is.null(useonly_t_MONITOR) && !is.na(useonly_t_MONITOR)) assert_that(is.character(useonly_t_MONITOR))
  # OData$dat.sVar[, c("g0.CAN.compare") := list(h_gN)] # should be identical to g0.CAN
  if (!is.null(type_intervened_TRT) && !is.na(type_intervened_TRT)) assert_that(is.character(type_intervened_TRT))
  if (!is.null(type_intervened_MONITOR) && !is.na(type_intervened_MONITOR)) assert_that(is.character(type_intervened_MONITOR))
  # print("CALLING IP WEIGHTS NOW"); print("intervened_TRT"); print(intervened_TRT)
  # ------------------------------------------------------------------------------------------
  # Probabilities of counterfactual interventions under observed (A,C,N) at each t
  # Combine the propensity score for observed (g0.C, g0.A, g0.N) with the propensity scores for interventions (gstar.C, gstar.A, gstar.N):
  # ------------------------------------------------------------------------------------------------------------------------------
  # (1) gstar.CENS: the indicator of not being censored.
  # (2) gstar.TRT: prob of following one treatment rule; and
  # (3) gstar.MONITOR prob following the monitoring regime; and
  # ------------------------------------------------------------------------------------------------------------------------------
  if (is.null(g_preds) || is.na(g_preds))
    stop("...cannot locate propensity scores in 'OData' object - must run fitPropensity(...) prior to calling this function")
  if (any(!(c("g0.A", "g0.C", "g0.N") %in% names(g_preds))))
    stop("... fatal error; propensity scores were not found in the input dataset, please re-run fitPropensity(...)")

  # indicator that the person is uncensored at each t (continuation of follow-up)
  gstar.CENS = as.integer(OData$eval_uncensored())
  # Likelihood P(A^*(t)=A(t)) under counterfactual intervention A^*(t) on A(t)
  gstar.TRT <- defineNodeGstarIPW(OData, intervened_TRT, nodes$Anodes, useonly_t_TRT, g_preds[["g0.A"]], OData$modelfit.gA, type_intervened_TRT)
  # Likelihood for monitoring P(N^*(t)=N(t)) under counterfactual intervention N^*(t) on N(t):
  gstar.MONITOR <- defineNodeGstarIPW(OData, intervened_MONITOR, nodes$Nnodes, useonly_t_MONITOR, g_preds[["g0.N"]], OData$modelfit.gN, type_intervened_MONITOR)

  # Save all likelihoods relating to propensity scores in separate dataset:
  # wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, nodes$Ynode, "g0.A", "g0.C", "g0.N", "g0.CAN"), with = FALSE]
  wts.DT <- OData$dat.sVar[, c(nodes$IDnode, nodes$tnode, nodes$Ynode), with = FALSE]
  wts.DT[,
    c("g0.A", "g0.C", "g0.N") := list(g_preds[["g0.A"]], g_preds[["g0.C"]], g_preds[["g0.N"]])][,
    c("g0.CAN") := g0.A * g0.C * g0.N]

  setkeyv(wts.DT, cols = c(nodes$IDnode, nodes$tnode))

  wts.DT[,
    "gstar.C" := gstar.CENS][,
    "gstar.A" := gstar.TRT][,
    "gstar.N" := gstar.MONITOR][,
    "gstar.CAN" := gstar.CENS * gstar.TRT * gstar.MONITOR] # Joint likelihoood for all 3 node types:

  ## Weights by time and cumulative weights by time:
  wts.DT[,"wt.by.t" := gstar.CAN / g0.CAN, by = eval(nodes$IDnode)]

  ## When ignore_tmin is specified set the wt.by.t to 1 for all t values that occur prior to time-points.
  ## NOTE: This is not the most efficient way to evaluate forward product for ignore_tmin, but it preserves the indexing
  ## of the original database, making it very easy to look-up correct weights for each observation row-index from the main dataset.
  if (!is.null(ignore_tmin)) wts.DT[eval(as.name(nodes$tnode)) < ignore_tmin, ("wt.by.t") := 1]
  if (!is.null(ignore_tmax)) wts.DT[eval(as.name(nodes$tnode)) > ignore_tmax, ("wt.by.t") := 1]
  if (reverse_wt_prod) {
    wts.DT[,"cum.IPAW" := rev(cumprod(rev(wt.by.t))), by = eval(nodes$IDnode)]
  } else {
    wts.DT[,"cum.IPAW" := cumprod(wt.by.t), by = eval(nodes$IDnode)]
  }

  if (trunc_weights < Inf) wts.DT[eval(as.name("cum.IPAW")) > trunc_weights, ("cum.IPAW") := trunc_weights]

  ## -------------------------------------------------------------------------------------------
  ## Calculate weight stabilization factor -- get emp P(followed rule at time t | followed rule up to now)
  ## -------------------------------------------------------------------------------------------
  if (eval_stabP) {
    nIDs <- OData$nuniqueIDs
    # THE ENUMERATOR: the total sum of subjects followed the rule gstar.A at t
    # THE DENOMINATOR: divide above by the total number of subjects who were still at risk of NOT FOLLOWING the rule at t
    # i.e., followed rule at t-1, assume at the first time-point EVERYONE was following the rule (so denominator = n)
    # (The total sum of all subjects who WERE AT RISK at t)
    # (FASTER) Version outside data.table, then merge back results:
    wts.DT[, "rule.follower.gCAN" := as.integer(cumprod(gstar.CAN) > 0), by = eval(nodes$IDnode)]
    n.follow.rule.t <- wts.DT[, list(N.follow.rule = sum(rule.follower.gCAN, na.rm = TRUE)), by = eval(nodes$tnode)]
    wts.DT[, "rule.follower.gCAN" := NULL]

    n.follow.rule.t[,
      N.risk := shift(N.follow.rule, fill = nIDs, type = "lag")][,
      stab.P := ifelse(N.risk > 0, N.follow.rule / N.risk, 1)][,
      cum.stab.P := cumprod(stab.P)]

    n.follow.rule.t[, c("N.risk", "stab.P") := list(NULL, NULL)]
    setkeyv(n.follow.rule.t, cols = nodes$tnode)

    wts.DT <- wts.DT[n.follow.rule.t, on = nodes$tnode]
    setkeyv(wts.DT, cols = c(nodes$IDnode, nodes$tnode))
  }

  ## remove person time observations with FUP above tmax
  if (!is.null(tmax)) wts.DT <- wts.DT[eval(as.name(nodes$tnode)) <= tmax, ]
  setkeyv(wts.DT, cols = c(nodes$IDnode, nodes$tnode))

  wts.DT[, "rule.name" := eval(as.character(rule_name))]

  attributes(wts.DT)[['getIPWeights_fun_call']] <- getIPWeights_fun_call
  attributes(wts.DT)[['intervened_TRT']] <- intervened_TRT
  attributes(wts.DT)[['intervened_MONITOR']] <- intervened_MONITOR
  return(wts.DT)
}

# ---------------------------------------------------------------------------------------
#' Direct (bounded) IPW estimator for the expected outcome over time (can be time-to-event or not).
#' @param wts_data \code{data.table} returned by a single call to \code{getIPWeights}.
#' Must contain the treatment/monitoring estimated IPTW weights for a SINGLE rule.
#' @param OData The object returned by function \code{fitPropensity}.
#' Contains the input data and the previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param weights (NOT IMPLEMENTED) Optional \code{data.table} with additional observation-time-specific weights.
#' Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param trunc_weights (NOT IMPLEMENTED) Specify the numeric weight truncation value.
#' All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @param return_wts Return the data.table with subject-specific IP weights as part of the output.
#' Note: for large datasets setting this to \code{TRUE} may lead to extremely large object sizes!
#' @return A data.table with bounded IPW estimates of risk and survival by time.
#' @example tests/examples/2_building_blocks_example.R
#' @export
directIPW <- function(wts_data, OData, weights, trunc_weights = 10^6, return_wts = FALSE) {
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode
  rule.name <- unique(wts_data[["rule.name"]])
  if (length(rule.name)>1)
    stop("wts_data must contain the weights for a single rule, found more than one unique rule name under in 'rule.name' column")

  ## EXTRACT RELEVANT INFORMATION
  ID.t.IPW.Y <- wts_data[,list(get(nodes$IDnode), get(t_name), cum.IPAW, get(Ynode))]
  names(ID.t.IPW.Y) <- c(nodes$IDnode, t_name, "cum.IPAW", Ynode)
  setkeyv(ID.t.IPW.Y, cols = c(nodes$IDnode, t_name))

  ## MAKE SURE EVERY PATIENT HAS AN ENTRY FOR EVERY TIME POINT BY LVCF
  ## version 1:
  # tmax <- wts_data[, max(get(t_name))]
  # tmax <- tmax - 1
  # UID <- wts_data[, unique(get(nodes$IDnode))]
  # all.ID.t <- as.data.table(cbind(rep(UID,each=(tmax+1)), rep(0:tmax,length(UID)) ))
  # names(all.ID.t) <- c(nodes$IDnode, t_name)
  # all.ID.t <- merge(all.ID.t, ID.t.IPW.Y, all.x=TRUE, by = c(nodes$IDnode, t_name))
  # all.ID.t[ , c("cum.IPAW", Ynode) := list(zoo::na.locf(cum.IPAW), zoo::na.locf(get(Ynode))), by = get(nodes$IDnode)]

  ## version 2:
  all.ID.t <- data.table::CJ(unique(wts_data[[nodes$IDnode]]), unique(wts_data[[t_name]]))
  colnames(all.ID.t) <- c(nodes$IDnode, t_name)
  setkeyv(all.ID.t, cols = c(nodes$IDnode, t_name))
  all.ID.t <- merge(all.ID.t, ID.t.IPW.Y, all.x=TRUE) # , by = c(nodes$IDnode, t_name)
  ## carry forward last known outcome / weights for survival outcomes (after time-to-event to maximum follow-up)
  if (any(is.na(all.ID.t[[Ynode]]))) all.ID.t[ , c("cum.IPAW", Ynode) := list(zoo::na.locf(cum.IPAW), zoo::na.locf(get(Ynode))), by = get(nodes$IDnode)]

  ## NUMERATOR OF BOUNDED IPW FOR SURVIVAL:
  numIPW <- all.ID.t[, {sum_Y_IPW = sum(get(Ynode)*cum.IPAW, na.rm = TRUE); list(sum_Y_IPW = sum_Y_IPW)}, by = eval(t_name)]
  ## DENOMINATOR OF BOUNDED IPW FOR SURVIVAL:
  denomIPW <- all.ID.t[, {sum_IPW = sum(cum.IPAW, na.rm = TRUE); list(sum_IPW = sum_IPW)}, by = eval(t_name)]
  ## BOUNDED IPW OF SURVIVAL (DIRECT):
  risk.t <- (numIPW[, "sum_Y_IPW", with = FALSE] / denomIPW[, "sum_IPW", with = FALSE])

  resultDT <- data.table(merge(numIPW, denomIPW, by = t_name))
  # resultDT <- data.table(est_name = "DirectBoundedIPW", merge(numIPW, denomIPW, by = t_name))
  resultDT[, ("St.directIPW") := (1 - risk.t[[1]])]
  resultDT[, "rule.name" := rule.name]

  ## STANDARDIZE THE NAME OF THE 'time' VARIABLE
  setnames(resultDT, t_name, "time")

  ## FILL IN THE GAP BETWEEN MINIMAL 'time' VALUE FROM OBSERVED DATA (IF THERE IS ANY) BY ADDING EXTRA ROWS IN THE BEGINING
  ## NOT NEEDED FOR DIRECT IPW -- THIS SHOULD BE ACCOMPLISHED AUTOMATICALLY
  # if (min(resultDT[["time"]]) > OData$min.t) {
  #   t_to_add <- (OData$min.t) : (min(resultDT[["time"]])-1)
  #   resultDT_addt <- merge(data.table(time = t_to_add), resultDT, by = "time", all.x = TRUE)[,
  #                             c("St.NPMSM", "St.KM") := 1][,
  #                             c("ht.NPMSM","ht.KM") := 0][,
  #                             "rule.name" := resultDT[["rule.name"]][1]]
  #   resultDT <- rbind(resultDT_addt, resultDT)
  #   setkeyv(resultDT, cols = "time")
  # }

  est_name <- "directIPW"
  resultDT <- cbind(est_name = est_name, resultDT)
  # resultDT <- data.frame(resultDT)

  attr(resultDT, "estimator_short") <- est_name
  attr(resultDT, "estimator_long") <- "Direct Bounded IPW"
  attr(resultDT, "trunc_weights") <- trunc_weights
  attr(resultDT, "rule_name") <- rule.name
  attr(resultDT, "time") <- resultDT[["time"]]

  result_object <- list(estimates = resultDT)

  if (return_wts) result_object[["wts_data"]] <- wts_data

  return(result_object)
}

## ---------------------------------------------------------------------------------------
## ***** TO DO *****: rename to survNPMSM_1rule(), internally called to estimate for single rule.
## ***** TO DO *****: if more then one rule exists in wts_data, then loop over each calling this function, combining results
## ---------------------------------------------------------------------------------------
#' Non-parametric (saturated) MSM for survival based on previously evaluated IP weights.
#' @param wts_data \code{data.table} returned by a single call to \code{getIPWeights}.
#' Must contain the treatment/monitoring estimated IPTW weights for a SINGLE rule.
#' @param OData The object returned by function \code{fitPropensity}.
#' Contains the input data and the previously fitted propensity score models for the exposure, monitoring and
#' right-censoring.
#' @param weights Optional \code{data.table} with additional observation-time-specific weights.
#' Must contain columns \code{ID}, \code{t} and \code{weight}.
#' The column named \code{weight} is merged back into the original data according to (\code{ID}, \code{t}).
#' @param trunc_weights Specify the numeric weight truncation value.
#' All final weights exceeding the value in \code{trunc_weights} will be truncated.
#' @param return_wts Return the data.table with subject-specific IP weights as part of the output.
#' Note: for large datasets setting this to \code{TRUE} may lead to extremely large object sizes!
#' @return A data.table with hazard and survival function estimates by time.
#' Also include the unadjusted Kaplan-Maier estimates.
#' @example tests/examples/2_building_blocks_example.R
#' @export
survNPMSM <- function(wts_data,
                      OData,
                      weights = NULL,
                      trunc_weights = 10^6,
                      return_wts = FALSE) {
  nodes <- OData$nodes
  t_name <- nodes$tnode
  Ynode <- nodes$Ynode

  # if (gvars$verbose) {print("trunc_weights"); print(trunc_weights)}

  rule.name <- unique(wts_data[["rule.name"]])
  if (length(rule.name)>1)
    stop("wts_data must contain the weights for a single rule, found more than one unique rule name under in 'rule.name' column")

  wts_data_used <- wts_data[,c(nodes$IDnode,nodes$tnode,nodes$Ynode,"cum.stab.P","cum.IPAW"), with = FALSE]
  setkeyv(wts_data_used, cols = c(nodes$IDnode, nodes$tnode))

  ## INITIALIZE WEIGHTED OUTCOME 'Wt.OUTCOME' TO Ynode:
  wts_data_used[, "Wt.OUTCOME" := get(nodes$Ynode)]

  ## ADD THE OBSERVATION-SPECIFIC WEIGHTS TO THE WEIGHTED OUTCOME, MERGE IN BY ID & t:
  wts_data_used <- process_opt_wts(wts_data_used, weights, nodes)

  ## ------------------------------------------------------------------------
  ## CRUDE HAZARD AND KM ESTIMATE OF SURVIVAL:
  ## ------------------------------------------------------------------------
  ht.crude <- wts_data_used[cum.IPAW > 0, {ht.KM = sum(Wt.OUTCOME, na.rm = TRUE) / .N; list(ht.KM = ht.KM)}, by = eval(t_name)][, St.KM := cumprod(1 - ht.KM)]
  setkeyv(ht.crude, cols = t_name)

  ## ------------------------------------------------------------------------
  ## IPW-ADJUSTED KM (SATURATED MSM):
  ## ------------------------------------------------------------------------
  ## MULTIPLY THE WEIGHT BY STABILIZATION FACTOR (NUMERATOR) (DOESN'T CHANGE THE ESTIMAND IN SATURATED MSMs, BUT MAKES WEIGHTS SMALLER):
  wts_data_used[, "cum.IPAW" := cum.stab.P * cum.IPAW]

  ## IF trunc_weights < Inf THEN TRUNCATE THE WEIGHTS:
  if (trunc_weights < Inf) wts_data_used[cum.IPAW > trunc_weights, cum.IPAW := trunc_weights]

  ## MULTIPLY THE OUTCOME BY CUMULATIVE WEIGHTS IN cum.IPAW:
  wts_data_used[, "Wt.OUTCOME" := Wt.OUTCOME * cum.IPAW]

  ## THE ENUMERATOR FOR THE HAZARD AT t: THE WEIGHTED SUM OF SUBJECTS WHO HAD EXPERIENCED THE EVENT AT t:
  sum_Ywt <- wts_data_used[, {sum_Y_IPW = sum(Wt.OUTCOME, na.rm = TRUE); list(sum_Y_IPW = sum_Y_IPW)}, by = eval(t_name)]
  setkeyv(sum_Ywt, cols = t_name)

  ## THE DENOMINATOR FOR THE HAZARD AT t: THE WEIGHTED SUM OF ALL SUBJECTS WHO WERE AT RISK AT t (EQUIVALENT TO SUMMING CUMULATIVE WEIGHTS cum.IPAW BY t):
  sum_Allwt <- wts_data_used[, {sum_all_IPAW = sum(cum.IPAW, na.rm = TRUE); list(sum_all_IPAW = sum_all_IPAW)}, by = eval(t_name)]
  setkeyv(sum_Allwt, cols = t_name)

  ## EVALUATE THE DISCRETE HAZARD ht AND SURVIVAL St OVER t
  St_ht_IPAW <- sum_Ywt[sum_Allwt][, ("ht.NPMSM") := sum_Y_IPW / sum_all_IPAW][, ("St.NPMSM") := cumprod(1 - ht.NPMSM)]
  St_ht_IPAW <- merge(St_ht_IPAW, ht.crude, all = TRUE)
  St_ht_IPAW[, "rule.name" := rule.name]

  ## STANDARDIZE THE NAME OF THE 'time' VARIABLE
  setnames(St_ht_IPAW, t_name, "time")

  ## FILL IN THE GAP BETWEEN MINIMAL 'time' VALUE FROM OBSERVED DATA (IF THERE IS ANY) BY ADDING EXTRA ROWS IN THE BEGINING
  if (min(St_ht_IPAW[["time"]]) > OData$min.t) {
    t_to_add <- (OData$min.t) : (min(St_ht_IPAW[["time"]])-1)
    St_ht_IPAW_addt <- merge(data.table(time = t_to_add), St_ht_IPAW, by = "time", all.x = TRUE)[,
                              c("St.NPMSM", "St.KM") := 1][,
                              c("ht.NPMSM","ht.KM") := 0][,
                              "rule.name" := St_ht_IPAW[["rule.name"]][1]]
    St_ht_IPAW <- rbind(St_ht_IPAW_addt, St_ht_IPAW)
    setkeyv(St_ht_IPAW, cols = "time")
  }

  est_name <- "NPMSM"
  resultDT <- cbind(est_name = est_name, St_ht_IPAW)

  attr(resultDT, "estimator_short") <- est_name
  attr(resultDT, "estimator_long") <- "NPMSM (Non-Parametric Marginal Structural Model) / AKME (IPW Adjusted Kaplan-Meier)"
  attr(resultDT, "trunc_weights") <- trunc_weights
  attr(resultDT, "rule_name") <- rule.name
  attr(resultDT, "time") <- resultDT[["time"]]
  # estimates <- resultDT
  # names(estimates) <- rule.name

  result_object <- list(estimates = resultDT)

  if (return_wts) result_object[["wts_data"]] <- wts_data_used

  return(result_object)
}
