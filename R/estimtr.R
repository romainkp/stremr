
#' @useDynLib estimtr
#' @import R6
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics axis barplot hist par text  legend plot
#' @importFrom methods is
#' @importFrom stats approx binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str
#' @importFrom stats as.formula glm na.exclude rbinom terms.formula
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

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}
# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

checkpkgs <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg %+% " package needed for this function to work. Please install it.", call. = FALSE)
    }
  }
}

# Bound g(A|W) probability within supplied bounds
bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds)
  x[x>max(bounds)] <- max(bounds)
  return(x)
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
                            "prediction from a rank-deficient fit may be misleading")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
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

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}

# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
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
    strat.C <- paste0(as.vector(get_outvars(regs)[1:Var_indx]) %+% " == " %+% gvars$noCENS.cat, collapse=" & ")
    curr_exprs <- get_subset_exprs(regs)[[Var_indx+1]]
    if (!is.null(curr_exprs)) {
      reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_expr = stringr::str_c(curr_exprs, " & ", strat.C))
    } else {
      reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_expr = strat.C)
    }
  }
  return(regs)
}
# Create subsetting expressions for a node (Anode, Cnode or Nnode)
# Named list with character expressions for subsetting. Each list item corresponds to one outcome in SingleRegressionFormClass
create_subset_expr <- function(outvars, stratify.EXPRS) {
  if (is.null(stratify.EXPRS)) {
    return(NULL)
  }
  Node_subset_expr <- vector(mode="list", length=length(outvars))
  names(Node_subset_expr) <- outvars
  assert_that(is.list(stratify.EXPRS))
  if (!all(outvars %in% names(stratify.EXPRS))) {
    stop("Could not locate the appropriate regression variable(s) within the supplied stratification list stratify.CENS, stratify.TRT or stratify.MONITOR." %+% "\n" %+%
          "The regression outcome variable(s) specified in gform.CENS, gform.TRT or gform.MONITOR were: ( '" %+% paste0(outvars, collapse=",") %+% "' )" %+% "\n" %+%
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
process_regforms <- function(regforms, default.reg, stratify.EXPRS = NULL, OData, sVar.map = NULL, factor.map = NULL, censoring = FALSE) {
  using.default <- FALSE
  if (missing(regforms)) {
    using.default <- TRUE
    regforms <- default.reg
  }
  if (!is.null(stratify.EXPRS)) assert_that(is.list(stratify.EXPRS))
  outvars <- vector(mode="list", length=length(regforms))
  predvars <- vector(mode="list", length=length(regforms))
  regs <- vector(mode="list", length=length(regforms))

  for (idx in seq_along(regforms)) {
    res <- process_regform(as.formula(regforms[[idx]]), sVar.map = sVar.map, factor.map = factor.map)
    outvars[[idx]] <- res$outvars
    if (!is.null(res$predvars)) predvars[[idx]] <- res$predvars
    names(outvars)[idx] <- names(predvars)[idx] <- paste0(outvars[[idx]], collapse="+")
    if (using.default && gvars$verbose)
      message("Using the default regression formula: " %+% paste0(outvars[[idx]], collapse="+") %+% " ~ " %+% paste0(predvars[[idx]], collapse="+"))
      # browser()
      subset_expr <- create_subset_expr(outvars = res$outvars, stratify.EXPRS = stratify.EXPRS)
      regobj <- SingleRegressionFormClass$new(outvar = res$outvars, predvars = res$predvars, outvar.class = OData$type.sVar[res$outvars],
                                              subset_vars = NULL, subset_exprs = subset_expr, censoring = censoring)
      regs[[idx]] <- regobj
  }
  class(regs) <- c(class(regs), "ListOfRegressionForms")
  if (censoring) regs <- stratify_by_uncensored(regs)
  return(list(outvars = outvars, predvars = predvars, regs = regs))
}

#---------------------------------------------------------------------------------
# MAIN estimtr FUNCTION
#---------------------------------------------------------------------------------
#' Estimate Survival with Interventions on Exposure and MONITORing Process in Right Censored Longitudinal Data.
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the treatment/exposure and monitoring process.
#'  Implements the \strong{IPW} (Inverse Probability-Weighted or Horvitz-Thompson) estimator of the discrete survival hazard function which is mapped into survival function.
#' @param data Observed data in long format. Should be a \code{data.frame} with named columns, containing the time-varying covariates (\code{covars}),
#'  the right-censoring event indicator(s) (\code{CENS}), the exposure variable(s) (\code{TRT}), the monitoring process variable(s) (\code{MONITOR})
#'  and the survival OUTCOME variable (\code{OUTCOME}).
# @param estimators (NOT IMPLEMENTED) Character vector with estimator names.
#' @param ID Unique subject identifier variable in the input data.
#' @param t The name of the time/period variable in the input data.
#' @param covars Time varying and baseline covariates. This argument does not need to be specified, by default all variables
#' that are not in \code{ID}, \code{t}, \code{CENS}, \code{TRT}, \code{MONITOR} and \code{OUTCOME} will be considered as covariates.
#' @param CENS CENSoring variable(s) in the input data.
#' Each separate variable specified in \code{CENS} can be either binary (0/1 valued integer) or categorical (integer).
#' For binary indicators of CENSoring, the value of 1 indicates the CENSoring or end of follow-up event (this cannot be changed).
#' For categorical CENSoring variables, by default the value of 0 indicates no CENSoring / continuation of follow-up and other
#' values indicate different reasons for CENSoring.
#' Use the argument \code{noCENS.cat} to change the reference (continuation of follow-up) category from default 0 to any other value.
#' (NOTE: Changing \code{noCENS.cat} has zero effect on coding of the binary CENSoring variables, those have to always use 1 to code the CENSoring event).
#' Note that factors are not allowed in \code{CENS}.
#' @param TRT Exposure/treatment variable(s) in input data.
#' @param MONITOR Monitoring variable(s) in input data.
#' @param OUTCOME  Survival OUTCOME variable name (column name in \code{data}).
#' @param noCENS.cat The level (integer) that indicates CONTINUATION OF FOLLOW-UP for ALL censoring variables. Defaults is 0.
#' Use this to modify the default reference category (no CENSoring / continuation of follow-up)
#' for variables specifed in \code{CENS}.
#' @param gform.TRT  Regression formula(s) for propensity score for the exposure/treatment(s): P(A(t) | W). See Details.
#' @param gform.CENS  Regression formula(s) for estimating the propensity score for the censoring mechanism: P(C(t) | W). See Details.
#' @param gform.MONITOR  Regression formula(s) for estimating the propensity score for the MONITORing process: P(N(t) | W). See Details.
# @param hform.g0 Regression formula for estimating the conditional density of P(\code{sA} | \code{sW}) under \code{g0}
#' (the observed exposure mechanism), When omitted the regression is defined by \code{sA~sW}, where \code{sA}
#  are all summary measures defined by argument \code{sA} and \code{sW} are all baseline summary measures defined by argument \code{sW}.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(estimtr.verbose=TRUE)}.
#' @param optPars A named list of additional optional parameters to be passed to \code{estimtr}, such as
#'  \code{alpha}, \code{lbound}, \code{family}, \code{YnodeDET},
#'  \code{h_g0_SummariesModel} and \code{h_gstar_SummariesModel}. See Details below for the description of each parameter.
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
#' Some of the parameters that control the estimation in \code{estimtr} can be set by calling the function \code{\link{estimtr_options}}.
#'
#' Additional parameters can be also specified as a named list \code{optPars} argument of the \code{estimtr} function.
#' The items that can be specified in \code{optPars} are:
#' \itemize{
#'
#' \item \code{alpha} - alpha-level for CI calculation (0.05 for 95% CIs);
#'
#' \item \code{lbound} - One value for symmetrical bounds on P(sW | sW).
#'
#' \item \code{family} - Family specification for regression models, defaults to binomial (CURRENTLY ONLY BINOMIAL
#'  FAMILY IS IMPLEMENTED).
#
#(REMOVED)\item \code{n_MCsims} - Number of Monte-Carlo simulations performed, each of sample size \code{nrow(data)},
#    for generating new exposures under \code{f_gstar1} or \code{f_gstar2} (if specified) or \code{f_g0} (if specified).
#    These newly generated exposures are utilized when fitting the conditional densities P(\code{sA}|\code{sW})
#    and when evaluating the substitution estimators \strong{GCOMP} and \strong{TMLE}
#    under stochastic interventions \code{f_gstar1} or \code{f_gstar2}.
#(REMOVED) \item \code{h_g0_SummariesModel} - Previously fitted model for P(\code{sA}|\code{sW}) under observed exposure mechanism \code{g0},
#    returned by the previous runs of the \code{estimtr} function.
#    This has to be an object of \code{SummariesModel} \pkg{R6} class. When this argument is specified, all predictions
#    P(\code{sA}=\code{sa}|\code{sW}=\code{sw}) under \code{g0} will be based on the model fits provided by this argument.
#(REMOVED) \item \code{h_gstar_SummariesModel} - Previously fitted model for P(\code{sA}|\code{sW}) under (stochastic) intervention
#    specified by \code{f_gstar1} or \code{f_gstar2}. Also an object of \code{SummariesModel} \pkg{R6} class.
#    When this argument is specified, the predictions P(\code{sA}=\code{sa}|\code{sW}=\code{sw})
#    under \code{f_gstar1} or \code{f_gstar2} will be based on the model fits provided by this argument.
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
#' @return A named list with 3 items containing the estimation results for:
#'  \itemize{
#'  \item \code{EY_gstar1} - estimates of the mean counterfactual OUTCOME under (stochastic) intervention function \code{f_gstar1} \eqn{(E_{g^*_1}[Y])}.
#'  \item \code{EY_gstar2} - estimates of the mean counterfactual OUTCOME under (stochastic) intervention function \code{f_gstar2} \eqn{(E_{g^*_2}[Y])},
#'    \code{NULL} if \code{f_gstar2} not specified.
#'  \item \code{ATE} - additive treatment effect (\eqn{E_{g^*_1}[Y]} - \eqn{E_{g^*_2}[Y]}) under interventions \code{f_gstar1}
#'    vs. in \code{f_gstar2}, \code{NULL} if \code{f_gstar2} not specified.
#' }
#'
#' Each list item above is itself a list containing the items:
#'  \itemize{
#'  \item \code{estimates} - various estimates of the target parameter (network population counterfactual mean under
#'    (stochastic) intervention).
#'  \item \code{vars} - the asymptotic variance estimates for \strong{IPTW}.
#'  \item \code{CIs} - CI estimates at \code{alpha} level for \strong{IPTW}.
#'  \item \code{other.vars} - Placeholder for future versions.
# \item \code{h_g0_SummariesModel} - The model fits for P(\code{sA}|\code{sW}) under observed exposure mechanism
#    \code{g0}. This is an object of \code{SummariesModel} \pkg{R6} class.
#  \item \code{h_gstar_SummariesModel} - The model fits for P(\code{sA}|\code{sW}) under intervention \code{f_gstar1}
#    or \code{f_gstar2}. This is an object of \code{SummariesModel} \pkg{R6} class.
#' }
#'
#' Currently implemented estimators are:
#'  \itemize{
#'  \item \code{iptw} - IPTW
#' }
#' @seealso \code{\link{estimtr-package}} for the general overview of the package,
#' @example tests/examples/1_estimtr_example.R
#' @export
# ------------------------------------------------------------------------------------------------------------------------------
# TO DO:
# ------------------------------------------------------------------------------------------------------------------------------
# **** TO DO: ****
# Write a helper method for SummaryModel class which goes down the nested tree of objects and obtains the model fits for each regression/outcome
# - NEED TO IMPLEMENT $get.fits() METHOD IN SummariesModel which recursively calls itself down the model tree until it reaches BinOutModel and returns its fit (regression + coefficients)
#   The method needs to appropriately format the output based on several model predictions (for stratified, categorical or continuous outcome)
# **** TO DO: ****
# IF A CERTAIN REGRESSION FORMULA / INTERVENTION NODE IS NOT SPECIFIED CREATE A DUMMY CLASS WHICH WOULD ALWAYS PUT MASS 1 ON THE OBERVED O
# TO DO: CONSIDER NOT THROWING AN ERROR WHEN stratify.VAR list is unnamed for cases where VAR is univariate (only one variable name)
# - Implement automatic function calling for gstar.TRT & gstar.MONITOR based on user-specified rule functions
#   If its a list of functions or if function returns more than one rule, apply the whole estimation procedure to each combination of TRT/MONITORING rules
# - Save the weights at each t and save the cummulative weights for all observations who were following the rule (g.CAN(O_i)>0)
# - Check that CENSor is either binary (integer or convert to integer) or categorical (integer or convert to integer)
# - look into g-force optimized functions for data.table: https://github.com/Rdatatable/data.table/issues/523
# ------------------------------------------------------------------------------------------------------------------------------
# - (IMPLEMENTED) Flip the CENSoring indicator for categorical CENS to make sure the reference category (noCENS.cat) IS ALWAYS CODED AS LAST
# - (IMPLEMENTED) When CENS[i] is binary and length(CENS)>1:
  # (1) Specify subset rule for i>1: (CENS[1]==0 & CENS[2]==0 & ... & CENS[i-1]==0)
  # (2) An alternative: collapse CENS into a categorical (automatically), based on the ordering in CENS. Then the fitting of categoricals will perform all the subsetting correctly.
  # (3) Alternative: set the indicators of missingness in the right way for CENS[i] if any CENS[1], ..., CENS[i-1] are 1.
# - (IMPLEMENTED) Stratification - allows K models on the SAME OUTCOME by stratifying rule
  # (1) User specified rule function creates strata. (stratify.CENS, stratify.TRT, stratify.MONITOR) Note that if the rule is based on data.table syntax it will be VERY FAST!
  # A nice trick would be to be able to AUTOMATICALLY convert logical subset expressions to data.table statements -> Its possible with some meta-programming and parsing
  # (2) These subsets (logical vectors) define K regressions, one regression model for each subset expression
  # Can specify K regressions in gform.CENS/gform.TRT/gform.MONITOR. If only one regression is specified it will be aplied to ALL stratas.
  # Otherwise stratas should be specified as a named list of K items
estimtr <- function(data, ID = "Subj_ID", t = "time_period",
                              covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                              gform.CENS, gform.TRT, gform.MONITOR, noCENS.cat = 0L,
                              stratify.CENS = NULL, stratify.TRT = NULL,
                              stratify.MONITOR = NULL,
                              gstar.TRT = NULL, gstar.MONITOR = NULL,
                              verbose = FALSE, optPars = list()) {

  gvars$verbose <- TRUE
  gvars$noCENS.cat <- noCENS.cat
  # if (verbose) {
    message("Running with the following setting: ");
    str(gvars$opts)
  # }

  gform.CENS.default <- "Cnodes ~ Lnodes"
  gform.TRT.default <- "Anodes ~ Lnodes"
  gform.MONITOR.default <- "Nnodes ~ Anodes + Lnodes"
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

  # ---------------------------------------------------------------------------
  # DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
  # Might expand this in the future to allow defining arbitrary summaries
  # ---------------------------------------------------------------------------
  lagnodes <- c(nodes$Cnodes, nodes$Anodes, nodes$Nnodes)
  newVarnames <- lagnodes %+% ".tminus1"

  print(str(lagnodes))
  print(OData$dat.sVar)

  OData$dat.sVar[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=get(nodes$ID), .SDcols=(lagnodes)]

  for (Cnode in nodes$Cnodes) CheckVarNameExists(OData$dat.sVar, Cnode)
  for (Anode in nodes$Anodes) CheckVarNameExists(OData$dat.sVar, Anode)
  for (Nnode in nodes$Nnodes) CheckVarNameExists(OData$dat.sVar, Nnode)
  for (Ynode in nodes$Ynode)  CheckVarNameExists(OData$dat.sVar, Ynode)
  for (Lnode in nodes$Lnodes) CheckVarNameExists(OData$dat.sVar, Lnode)

  # ------------------------------------------------------------------------------------------------
  # Process the input formulas and stratification settings;
  # Define regression classes for g.C, g.A, g.N and put them in a single list of regressions.
  # ------------------------------------------------------------------------------------------------
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
  summeas.g0 <- newsummarymodel(reg = ALL_g_regs, DatNet.sWsA.g0 = OData)
  summeas.g0$fit(data = OData)
  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...))
  h_gN <- summeas.g0$predictAeqa(newdata = OData)

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
  # Observed likelihood of (A,C,N) at each t
  # ------------------------------------------------------------------------------------------
  summeas.gC <- summeas.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gC")]]
  # summeas.gC$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$getPsAsW.models()
  summeas.gA <- summeas.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gA")]]
  # summeas.gA$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$getPsAsW.models()
  summeas.gN <- summeas.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gN")]]
  # summeas.gN$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$getfit

  g0.A <- summeas.gA$getcumprodAeqa()
  g0.C <- summeas.gC$getcumprodAeqa()
  g0.N <- summeas.gN$getcumprodAeqa()
  N_IDs <- length(unique(OData$dat.sVar[[ID]])) # Total number of observations
  OData$dat.sVar[, c("g0.A", "g0.C", "g0.N", "g0.CAN") := list(g0.A, g0.C, g0.N, g0.A*g0.C*g0.N)]
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
  OData$dat.sVar[, "gstar.C" := as.integer(rowSums(.SD) == eval(noCENS.cat)), .SDcols = CENS]

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
  OData$dat.sVar[, "wt.by.t" := gstar.CAN / g0.CAN, by = eval(ID)][, "cumm.IPAW" := cumprod(wt.by.t), by = eval(ID)]
  # OData$dat.sVar[1:100,]

  # -------------------------------------------------------------------------------------------
  # Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
  # NOTE: DO THIS AT THE VERY BEGINNING INSTEAD????
  # DEPENDS ON STRUCTURE OF THE DATA AND IF ANY C events ACTUALLY OCCURRED IN LAST ROW -> THIS IS THE CASE FOR ADMINISTRATIVE CENSORING
  # -------------------------------------------------------------------------------------------
  shifted.OUTCOME <- OUTCOME%+%".tplus1"
  OData$dat.sVar[, (shifted.OUTCOME) := shift(get(OUTCOME), n = 1L, type = "lead"), by = eval(ID)]
  OData$dat.sVar <- OData$dat.sVar[!is.na(get(shifted.OUTCOME)), ] # drop and over-write previous data.table, removing last rows.

  # multiply the shifted outcomes by the current (cummulative) weight cumm.IPAW:
  OData$dat.sVar[, "Wt.OUTCOME" := get(shifted.OUTCOME)*cumm.IPAW]
  # OData$dat.sVar[101:200, ]
  # Row indices for all subjects at t who had the event at t+1 (NOT USING)
  # row_idx_outcome <- OData$dat.sVar[, .I[get(shifted.OUTCOME) %in% 1L], by = eval(ID)][["V1"]]

  # THE ENUMERATOR FOR THE HAZARD AT t: the weighted sum of subjects who had experienced the event at t:
  sum_Ywt <- OData$dat.sVar[, .(sum_Y_IPAW=sum(Wt.OUTCOME)), by = eval(t)]; setkeyv(sum_Ywt, cols=t)
  # THE DENOMINATOR FOR THE HAZARD AT t: The weighted sum of all subjects who WERE AT RISK at t:
  # (equivalent to summing cummulative weights cumm.IPAW by t)
  sum_Allwt <- OData$dat.sVar[, .(sum_all_IPAW=sum(cumm.IPAW)), by = eval(t)]; setkeyv(sum_Allwt, cols=t)
  # EVALUATE THE DISCRETE HAZARD ht AND SURVIVAL St OVER t
  St_ht_IPAW <- sum_Ywt[sum_Allwt][, "ht" := sum_Y_IPAW / sum_all_IPAW][, c("m1ht", "St") := .(1-ht, cumprod(1-ht))]

return(list(IPW_estimates = data.frame(St_ht_IPAW), dataDT = OData$dat.sVar, model.R6.fits = summeas.g0, data.R6.object = OData))
}












