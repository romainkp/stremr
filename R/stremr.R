#' @useDynLib stremr
#' @import R6
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics axis barplot hist par text  legend plot
#' @importFrom methods is
#' @importFrom stats approx binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str capture.output
#' @importFrom stats as.formula glm na.exclude rbinom terms.formula pnorm quasibinomial
NULL

#-----------------------------------------------------------------------------
# Class Membership Tests
#-----------------------------------------------------------------------------
is.DataStorageClass <- function(DataStorageClass) "DataStorageClass"%in%class(DataStorageClass)

#-----------------------------------------------------------------------------
# General utilities / Global Vars
#-----------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
is.integerish <- function (x) is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))

# # Return the left hand side variable of formula f as a character
# LhsVars <- function(f) {
#   f <- as.formula(f)
#   return(as.character(f[[2]]))
# }
# # Return the right hand side variables of formula f as a character vector
# RhsVars <- function(f) {
#   f <- as.formula(f)
#   return(all.vars(f[[3]]))
# }
# # Bound g(A|W) probability within supplied bounds
# bound <- function(x, bounds){
#   x[x<min(bounds)] <- min(bounds)
#   x[x>max(bounds)] <- max(bounds)
#   return(x)
# }

checkpkgs <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg %+% " package needed for this function to work. Please install it.", call. = FALSE)
    }
  }
}

# returns NULL if no factors exist, otherwise return the name of the factor variable(s)
CheckExistFactors <- function(data) {
  testvec <- unlist(lapply(data, is.factor))
  if (any(testvec)) {
    return(names(data)[which(testvec)])
  } else {
    return(NULL)
  }
}

# throw exception if 1) varname doesn't exist; 2) more than one varname is matched
CheckVarNameExists <- function(data, varname) {
  idvar <- names(data) %in% varname
  if (sum(idvar) < 1) stop("variable name " %+% varname %+% " not found in data input")
  if (sum(idvar) > 1) stop("more than one column in the input data has been matched to name "
                            %+% varname %+% ". Consider renaming some of the columns: " %+%
                            paste0(names(data)[idvar], collapse=","))
  return(invisible(NULL))
}

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                            "prediction from a rank-deficient fit may be misleading",
                            "non-integer #successes in a binomial glm!",
                            "the matrix is either rank-deficient or indefinite")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

#---------------------------------------------------------------------------------
# SPECIFYING regressions for gC, gA & gN
#---------------------------------------------------------------------------------
get_vars_fromlist <- function(varname, sVar.map) {
  if (varname %in% names(sVar.map)) {
    as.vector(sVar.map[[varname]])
  } else {
    varname
  }
}
# Parse the formulas for summary measure names and create a map to actual covariate names in sA & sW
process_regform <- function(regform, sVar.map = NULL, factor.map = NULL) {
  # regform1 <- as.formula("N ~ 1")
  # Getting predictors (sW names):
  regformterms <- terms(regform)
  sW.names <- attributes(regformterms)$term.labels
  sW.names.alt <- colnames(attributes(regformterms)$factors)
  assert_that(all(sW.names == sW.names.alt))
  # Getting OUTCOMEs (sA names):
  (out.var <- deparse(attributes(regformterms)$variables[[2]])) # LHS character string

  out.vars.form <- as.formula(". ~ " %+% out.var)
  out.vars.terms <- terms(out.vars.form)
  sA.names <- attributes(out.vars.terms)$term.labels

  outvars <- unlist(lapply(sA.names, get_vars_fromlist, sVar.map))
  predvars <- unlist(lapply(sW.names, get_vars_fromlist, sVar.map))
  # in case some factors were also involved (these will then be replaced only on a second iteration)
  predvars <- unlist(lapply(predvars, get_vars_fromlist, factor.map))
  return(list(outvars = outvars, predvars = predvars))
}

# Loop through a list of SingleRegressionFormClass objects and their outvars as if it was one long list of outvars and create the subsetting expressions
# This uses S3 method dispatch on object ListOfRegressionForms
stratify_by_uncensored <- function(regs) {
  for (Var_indx in seq_along(get_outvars(regs)[-1])) {
    curr_outvar <- get_outvars(regs)[Var_indx+1]
    curr_exprs <- get_subset_exprs(regs)[[Var_indx+1]]
    prev_outvars <- as.vector(unique(get_outvars(regs)[1:Var_indx]))
    prev_outvars_to_condition <- prev_outvars[!(prev_outvars %in% curr_outvar)]
    if (length(prev_outvars_to_condition) > 0) {
      strat.C <- paste0(prev_outvars_to_condition %+% " == " %+% gvars$noCENScat, collapse=" & ")
      if (!is.null(curr_exprs)) {
        reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_exprs = stringr::str_c(curr_exprs, " & ", strat.C))
      } else {
        reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_exprs = strat.C)
      }
    }
  }
  return(regs)
}

# Create subsetting expressions for a node (Anode, Cnode or Nnode)
# Named list with character expressions for subsetting.
# Each list item corresponds to one outcome in SingleRegressionFormClass
create_subset_expr <- function(outvars, stratify.EXPRS) {
  if (is.null(stratify.EXPRS)) {
    return(NULL)
  }
  Node_subset_expr <- vector(mode="list", length=length(outvars))
  names(Node_subset_expr) <- outvars
  assert_that(is.list(stratify.EXPRS))
  if (!all(outvars %in% names(stratify.EXPRS))) {
    stop("Could not locate the appropriate regression variable(s) within the supplied stratification list stratify_CENS, stratify_TRT or stratify_MONITOR." %+% "\n" %+%
          "The regression outcome variable(s) specified in gform_CENS, gform_TRT or gform_MONITOR were: ( '" %+% paste0(outvars, collapse=",") %+% "' )" %+% "\n" %+%
          "However, the item names in the matching stratification list were: ( '" %+% paste0(names(stratify.EXPRS), collapse=",") %+% "' )"
          )
  }
  for (idx in seq_along(Node_subset_expr))
    if (!is.null(stratify.EXPRS[[outvars[idx]]]))
      Node_subset_expr[[idx]] <- stratify.EXPRS[[outvars[idx]]]
  return(Node_subset_expr)
}

# When several reg forms are specified (multivariate Anodes), process outvars into one vector and process predvars in a named list of vectors
# stratify.EXPRS - Must be a named list. One item (characeter vectors) per one outcome in regforms.
process_regforms <- function(regforms, default.reg, stratify.EXPRS = NULL, model_contrl = NULL, OData, sVar.map = NULL, factor.map = NULL, censoring = FALSE, outvar.class = NULL) {
  using.default <- FALSE
  if (missing(regforms)) {
    using.default <- TRUE
    regforms <- default.reg
  }
  if (!is.null(stratify.EXPRS)) assert_that(is.list(stratify.EXPRS))
  regs <- vector(mode="list", length=length(regforms))
  for (idx in seq_along(regforms)) {
    res <- process_regform(as.formula(regforms[[idx]]), sVar.map = sVar.map, factor.map = factor.map)
    if (using.default && gvars$verbose)
      message("Using the default regression formula: " %+% paste0(res$outvars, collapse="+") %+% " ~ " %+% paste0(res$predvars, collapse="+"))

      if (!is.null(outvar.class)) {
        outvar.class <- as.list(rep.int(outvar.class, length(res$outvars)))
        names(outvar.class) <- res$outvars
      } else {
        outvar.class <- OData$type.sVar[res$outvars]
        names(outvar.class) <- res$outvars
      }
      subset_exprs <- create_subset_expr(outvars = res$outvars, stratify.EXPRS = stratify.EXPRS)

      regobj <- RegressionClass$new(outvar = res$outvars, predvars = res$predvars, outvar.class = outvar.class,
                                    subset_vars = NULL, subset_exprs = subset_exprs, model_contrl = model_contrl,
                                    censoring = censoring)
      regs[[idx]] <- regobj
      outvar.class <- NULL
  }
  class(regs) <- c(class(regs), "ListOfRegressionForms")
  if (censoring) regs <- stratify_by_uncensored(regs)
  return(regs)
}

#---------------------------------------------------------------------------------
# MAIN stremr FUNCTION
#---------------------------------------------------------------------------------
#' Estimate Survival with Interventions on Exposure and MONITORing Process in Right Censored Longitudinal Data.
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the treatment/exposure and monitoring process.
#'  Implements the \strong{IPW} (Inverse Probability-Weighted or Horvitz-Thompson) estimator of the discrete survival hazard function which is mapped into survival function.
#' @param data Input data in long format. Can be a \code{data.frame} or a \code{data.table} with named columns, containing the time-varying covariates (\code{covars}),
#'  the right-censoring event indicator(s) (\code{CENS}), the exposure variable(s) (\code{TRT}), the monitoring process variable(s) (\code{MONITOR})
#'  and the survival OUTCOME variable (\code{OUTCOME}).
# @param estimators (NOT IMPLEMENTED) Character vector with estimator names.
#' @param ID Unique subject identifier column name in \code{data}.
#' @param t.name The name of the time/period variable in \code{data}.
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
#' @param OUTCOME  A column name in \code{data} for the survival OUTCOME variable name, code as 1 for the outcome event.
#' @param noCENScat The level (integer) that indicates CONTINUATION OF FOLLOW-UP for ALL censoring variables. Defaults is 0.
#' Use this to modify the default reference category (no CENSoring / continuation of follow-up)
#' for variables specifed in \code{CENS}.
#' @param gform_TRT  Regression formula(s) for propensity score for the exposure/treatment(s): P(A(t) | W). See Details.
#' @param gform_CENS  Regression formula(s) for estimating the propensity score for the censoring mechanism: P(C(t) | W). See Details.
#' @param gform_MONITOR  Regression formula(s) for estimating the propensity score for the MONITORing process: P(N(t) | W). See Details.
# @param hform.g0 Regression formula for estimating the conditional density of P(\code{sA} | \code{sW}) under \code{g0}
#' (the observed exposure mechanism), When omitted the regression is defined by \code{sA~sW}, where \code{sA}
#  are all summary measures defined by argument \code{sA} and \code{sW} are all baseline summary measures defined by argument \code{sW}.
#' @param stratify_CENS A named list with one item per variable in \code{CENS}.
#' Each list item is a character vector of stratification subsets for the corresponding variable in \code{CENS}.
#' @param stratify_TRT A named list with one item per variable in \code{TRT}.
#' Each list item is a character vector of stratification subsets for the corresponding variable in \code{TRT}.
#' @param stratify_MONITOR A named list with one item per variable in \code{MONITOR}.
#' Each list item is a character vector of stratification subsets for the corresponding variable in \code{MONITOR}.
#' @param intervened_TRT Column name in \code{data} containing the counterfactual probabilities of following a specific treatment regimen.
#' @param intervened_MONITOR Column name in \code{data} containing the counterfactual probabilities of following a specific monitoring regimen.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#' @param optPars A named list of additional optional parameters to be passed to \code{stremr}, such as
#'  \code{alpha}, \code{lbound}, \code{family}, \code{YnodeDET},
#'  \code{h_g0_GenericModel} and \code{h_gstar_GenericModel}. See Details below for the description of each parameter.
# (REMOVED) \code{n_MCsims}
#((NOT IMPLEMENTED)) @param Q.SL.library SuperLearner libraries for OUTCOME, Q
#((NOT IMPLEMENTED)) @param g.SL.library SuperLearner libraries for treatment mechanism, g
# @param sW Summary measures constructed from baseline covariates alone. This must be an object of class
#  \code{DefineSummariesClass} that is returned by calling the function \code{\link{def_sW}}.
# @param sA Summary measures constructed from exposures \code{Anodes} and baseline covariates. This must be an object of class
#  \code{DefineSummariesClass} that is returned by calling the function \code{\link{def_sW}}.
# @param Anodes Exposure (treatment) variable name (column name in \code{data}); exposures can be either binary, categorical or continuous.
#  This variable can be instead specified with argument \code{sA} by adding a call \code{+def_sA(Anodes="ExposureVarName")} to \code{sA}.
# @param AnodeDET Optional column name for indicators of deterministic values of exposures in \code{Anodes},
#  should be coded as (\code{TRUE}/\code{FALSE}) or (\code{1}/\code{0});
#  observations with \code{AnodeDET}=\code{TRUE}/\code{1} are assumed to have deterministically assigned exposures
# @param intervene1.sA
# @param f_gstar1 Either a function or a vector of counterfactual exposures. If a function, must return
#  a vector of counterfactual exposures evaluated based on the summary measures matrix (\code{sW,sA}) passed as a named
#  argument \code{"data"}, therefore, the function in \code{f_gstar1} must have a named argument \code{"data"} in its signature.
#  The interventions defined by \code{f_gstar1} can be static, dynamic or stochastic. If \code{f_gstar1} is specified as a
#  vector, it must be of length \code{nrow(data)} or 1 (constant treatment assigned to all observations).
#  See Details below and Examples in "EQUIVALENT WAYS OF SPECIFYING INTERVENTION \code{f_gstar1}" for demonstration.
# @param intervene2.sA
# @param f_gstar2 Either a function or a vector of counterfactual exposure assignments.
#  Used for estimating contrasts (average treatment effect) for two interventions, if omitted, only the average
#  counterfactual OUTCOME under intervention \code{f_gstar1} is estimated. The requirements for \code{f_gstar2}
#  are identical to those for \code{f_gstar1}.
#'
#' @section Details:
#'
#' The regression formalas in \code{Qform}, \code{hform.g0} and \code{hform.gstar} can include any summary measures names defined in
#'  \code{sW} and \code{sA}, referenced by their individual variable names or by their aggregate summary measure names.
#'  For example, \code{hform.g0 = "netA ~ netW"} is equivalent to
#'  \code{hform.g0 = "A + A_netF1 + A_netF2 ~ W + W_netF1 + W_netF2"} for \code{sW,sA} summary measures defined by
#'  \code{def_sW(netW=W[[0:2]])} and \code{def_sA(netA=A[[0:2]])}.
#'
#' @section Additional parameters:
#'
#' Some of the parameters that control the estimation in \code{stremr} can be set by calling the function \code{\link{set_all_stremr_options}}.
#'
#' Additional parameters can be also specified as a named list \code{optPars} argument of the \code{stremr} function.
#' The items that can be specified in \code{optPars} are:
#' \itemize{
#'
#' \item \code{alpha} - alpha-level for CI calculation (0.05 for 95% CIs);
#'
#' \item \code{lbound} - One value for symmetrical bounds on P(sW | sW).
#'
#' \item \code{family} - Family specification for regression models, defaults to binomial (CURRENTLY ONLY BINOMIAL
#'  FAMILY IS IMPLEMENTED).
#' }
#'
#' @section Specifying the counterfactual intervention function (\code{f_gstar1} and \code{optPars$f_gstar2}):
#'
#' The functions \code{f_gstar1} and \code{f_gstar2} can only depend on variables specified by the combined matrix
#'  of summary measures (\code{sW},\code{sA}), which is passed using the argument \code{data}. The functions should
#'  return a vector of length \code{nrow(data)} of counterfactual treatments for observations in the input data.
#'
#' @section IPTW estimator:
#' **********************************************************************
#'
#' \itemize{
#' \item As described in the following section, the first step is to construct an estimator \eqn{P_{g_N}(A(t) | L(t))}
#'    for the probability of exposure \eqn{P_{g_0}(A(t) | W(t))}.
#'
#' \item Based on the user specified stochastic intervention, we can also obtain \eqn{P_{g^*_N}(A^*(t) | L(t) }
#'
#' \item Combining the two probabilities forms the basis of the IPTW estimator,
#'    which is evaluated at the observed N data points \eqn{O_i=((L_i(t), A_i(t): t=0,...,K), Y_i), i=1,...,N} and is given by
#'    \deqn{\psi^{IPTW}_n = \sum_{i=1,...,N}{Y_i \frac{P_{g^*_N}(A^*(t)=A_i(t) | L(t)=L_i(t))}{P_{g_N}(A(t)=A_i(t) | L(t)=L_i(t))}}.}
#' }
#'
#' @return ...
#' @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/1_stremr_example.R
#' @export
# ------------------------------------------------------------------------------------------------------------------------------
# **** BUILDING BLOCKS ****
# ------------------------------------------------------------------------------------------------------------------------------
# - BLOCK 1: Process inputs and define OData R6 object
# - BLOCK 2: define regression models, define a single RegressionClass & fit the propensity score for observed data, summary.g0 g0 (C,A,N)
# - BLOCK 3: evaluate weights based gstar_TRT, gstar_MONITOR and observed propensity scores g0, the input is modelfits.g0 and OData object
# - BLOCK 4A: Non-parametric MSM for survival, no weight stabilization, input either single weights dataset or a list of weights datasets
# - BLOCK 4B: Saturated MSM pooling many regimens, includes weight stabilization and using closed-form soluaton for the MSM (can only do saturated MSM)
# - BLOCK 4C: Parametric MSM pooling many regimens, includes weight stabilization and parametric MSM (can include saturated MSM)
# - BLOCK 5: Builds a report with weight distributions, survival estimates, etc.
# ------------------------------------------------------------------------------------------------------------------------------
# **** TO DO: ****
# ------------------------------------------------------------------------------------------------------------------------------
# - CHANGE TO stremrOptions(); use it to set options(stremr.... ), which then is read internally by gvars at startup:
#   strOptions(strict.width = "no", digits.d = 3, vec.len = 4,
#            formatNum = function(x, ...)
#                        format(x, trim = TRUE, drop0trailing = TRUE, ...))
# - Allow looping over regimens to return regimen-specific non-zero weight datasets or list of such dataset (data.tables) that can be then all stacked and used for one MSM
# - Implement automatic function calling for gstar_TRT & gstar_MONITOR if its a function or a list of functions
# - When node name is "NULL" (not specified), do not fit a model for it. create a dummy class which would always put mass 1 on the oberved o
#   The method needs to appropriately format the output based on several model predictions (for stratified, categorical or continuous outcome)
# - Allow specification of counterfactual trt & monitor vaules / counterfactual probabilities of trt & monitor = 1. map automatically into rule follors/non-followers
# - Consider not throwing an error when stratify.VAR list is unnamed for cases where VAR is univariate (only one variable name)
#   If its a list of functions or if function returns more than one rule, apply the whole estimation procedure to each combination of TRT/MONITORING rules
# - Check that CENSor is either binary (integer or convert to integer) or categorical (integer or convert to integer)
# - look into g-force optimized functions for data.table: https://github.com/Rdatatable/data.table/issues/523
# ------------------------------------------------------------------------------------------------------------------------------
stremr <- function(data, ID = "Subj_ID", t.name = "time_period",
                  covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                  gform_CENS, gform_TRT, gform_MONITOR,
                  stratify_CENS = NULL, stratify_TRT = NULL, stratify_MONITOR = NULL,
                  intervened_TRT = NULL, intervened_MONITOR = NULL, noCENScat = 0L,
                  verbose = getOption("stremr.verbose"), optPars = list()) {
  # ------------------------------------------------------------------
  # - BLOCK 1: Process inputs and define OData R6 object
  # ------------------------------------------------------------------
  OData <- importData(data, ID, t.name, covars, CENS, TRT, MONITOR, OUTCOME, noCENScat, verbose)
  # ------------------------------------------------------------------
  # - BLOCK 2: define regression models, define a single RegressionClass & fit the propensity score for observed data, summary.g0 g0 (C,A,N)
  # ------------------------------------------------------------------
  OData <- fitPropensity(OData, gform_CENS, gform_TRT, gform_MONITOR, stratify_CENS, stratify_TRT, stratify_MONITOR)
  # ---------------------------------------------------------------------------------------
  # - BLOCK 3: evaluate weights based gstar_TRT, gstar_MONITOR and observed propensity scores g0, the input is modelfits.g0 and OData object
  # ---------------------------------------------------------------------------------------
  wts.DT <- getIPWeights(OData, intervened_TRT, intervened_MONITOR)
  # ---------------------------------------------------------------------------------------
  # - BLOCK 4A: Non-parametric MSM for survival, with weight stabilization, input either single weights dataset or a list of weights datasets,
  # Each dataset containing weights non-zero weights for single regimen
  # ---------------------------------------------------------------------------------------
  IPW_estimates <- survNPMSM(wts.DT, OData)
  # ---------------------------------------------------------------------------------------
  # - BLOCK 5: Builds a report with weight distributions, survival estimates, etc.
  # ---------------------------------------------------------------------------------------
  # .... TO DO ....

return(list(IPW_estimates = IPW_estimates$IPW_estimates, wts_data = IPW_estimates$wts_data, dataDT = OData$dat.sVar, modelfits.g0.R6 = OData$modelfits.g0, OData.R6 = OData))
}
