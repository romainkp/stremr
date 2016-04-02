
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
# SPECIFYING regressions for g.C, g.A & g.N
#---------------------------------------------------------------------------------
get_vars_fromlist <- function(varname, sVar.map) {
  if (varname %in% names(sVar.map)) {
    as.vector(sVar.map[[varname]])
  } else {
    varname
  }
}
# Parse the formulas for summary measure names and create a map to actual covariate names in sA & sW
process_regform <- function(regform, sVar.map = NULL) {
  # Getting predictors (sW names):
  regformterms <- terms(regform)
  sW.names <- attributes(regformterms)$term.labels
  sW.names.alt <- colnames(attributes(regformterms)$factors)
  assert_that(all(sW.names == sW.names.alt))
  # Getting OUTCOMEs (sA names):
  out.var <- rownames(attributes(regformterms)$factors)[1] # character string
  out.vars.form <- as.formula(". ~ " %+% out.var)
  out.vars.terms <- terms(out.vars.form)
  sA.names <- attributes(out.vars.terms)$term.labels

  outvars <- unlist(lapply(sA.names, get_vars_fromlist, sVar.map))
  predvars <- unlist(lapply(sW.names, get_vars_fromlist, sVar.map))
  return(list(outvars = outvars, predvars = predvars))
}

# When several reg forms are specified (multivariate Anodes), process outvars into one vector and process predvars in a named list of vectors
process_regforms <- function(regforms, default.reg, sVar.map = NULL) {
  using.default <- FALSE
  if (missing(regforms)) {
    using.default <- TRUE
    regforms <- default.reg
  }

  outvars <- vector(mode="list", length=length(regforms))
  predvars <- vector(mode="list", length=length(regforms))

  for (idx in seq_along(regforms)) {
    res <- process_regform(as.formula(regforms[[idx]]), sVar.map = sVar.map)
    outvars[[idx]] <- res$outvars
    predvars[[idx]] <- res$predvars
    names(outvars)[idx] <- names(predvars)[idx] <- paste0(outvars[[idx]], collapse="+")
    if (using.default && gvars$verbose)
      message("Using the default regression formula: " %+% paste0(outvars[[idx]], collapse="+") %+% " ~ " %+% paste0(predvars[[idx]], collapse="+"))
  }
  return(list(outvars = outvars, predvars = predvars))
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
estimtr <- function(data, ID = "Subj_ID", t = "time_period",
                              covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                              gform.CENS, gform.TRT, gform.MONITOR, noCENS.cat = 0L,
                              stratify.CENS = NULL, stratify.TRT = NULL, stratify.MONITOR = NULL, verbose = FALSE, optPars = list()) {

  gvars$verbose <- TRUE
  # if (verbose) {
    message("Running with the following setting: ");
    str(gvars$opts)
    # message("Running tmlenet with the following settings from optPars arg of tmlenet(): ");
    # str(optPars)
  # }

  gform.CENS.default <- "Cnodes ~ Lnodes"
  gform.TRT.default <- "Anodes ~ Lnodes"
  gform.MONITOR.default <- "Nnodes ~ Anodes + Lnodes"
  if (missing(covars)) { # define time-varing covars (L) as everything else in data besides these vars
    covars <- setdiff(colnames(data), c(ID, t, CENS, TRT, MONITOR, OUTCOME))
  }
  # The ordering of variables in this list is the assumed temporal order!
  nodes <- list(Lnodes = covars, Cnodes = CENS, Anodes = TRT, Nnodes = MONITOR, Ynode = OUTCOME)
  OData <- DataStorageClass$new(Odata = data, nodes = nodes, noCENS.cat = noCENS.cat)

  for (Cnode in nodes$Cnodes) CheckVarNameExists(OData$dat.sVar, Cnode)
  for (Anode in nodes$Anodes) CheckVarNameExists(OData$dat.sVar, Anode)
  for (Nnode in nodes$Nnodes) CheckVarNameExists(OData$dat.sVar, Nnode)
  for (Ynode in nodes$Ynode)   CheckVarNameExists(OData$dat.sVar, Ynode)
  for (Lnode in nodes$Lnodes)  CheckVarNameExists(OData$dat.sVar, Lnode)

  g.C.sVars <- process_regforms(regforms = gform.CENS, default.reg = gform.CENS.default, sVar.map = nodes)
  g.A.sVars <- process_regforms(regforms = gform.TRT, default.reg = gform.TRT.default, sVar.map = nodes)
  g.N.sVars <- process_regforms(regforms = gform.MONITOR, default.reg = gform.MONITOR.default, sVar.map = nodes)

  # put all three regression models (C,A,N) into one regression object!
  all_outVar_nms <- c(g.C.sVars$outvars, g.A.sVars$outvars, g.N.sVars$outvars)
  all_predVar_nms <- c(g.C.sVars$predvars, g.A.sVars$predvars, g.N.sVars$predvars)
  all_subsets_expr <- lapply(all_outVar_nms, function(var) lapply(var, function(var) {var}))
  all_outVar_class <- lapply(all_outVar_nms, function(outVar_nm) OData$type.sVar[outVar_nm])

  ALL.g.regs <- RegressionClass$new(sep_predvars_sets = TRUE,
                                    outvar.class = all_outVar_class,
                                    outvar = all_outVar_nms,
                                    predvars = all_predVar_nms,
                                    subset = all_subsets_expr)
  ALL.g.regs$S3class <- "generic"
  # using S3 method dispatch on ALL.g.regs:
  summeas.g0 <- newsummarymodel(reg = ALL.g.regs, DatNet.sWsA.g0 = OData)

  summeas.g0$fit(data = OData)
  h_gN <- summeas.g0$predictAeqa(newdata = OData)

  # browser()
  # ALL.g.regs <- RegressionClass$new(outvar = nodes$Ynode,
  #                                                         predvars = Q.sVars$predvars[[1]],
  #                                                         subset = !determ.Q,
  #                                                         ReplMisVal0 = TRUE)
  # m.Q.init <- BinOutModel$new(glm = FALSE, reg = Qreg)$fit(data = OData)

  # - Check that CENSor is either binary (integer or convert to integer) or categorical (integer or convert to integer)

  # - Flip the CENSoring indicator for categorical CENS to make sure the reference category (noCENS.cat) IS ALWAYS CODED AS LAST

  # - When CENS[i] is binary and length(CENS)>1:
    # (1) Specify subset rule for i>1: (CENS[1]==0 & CENS[2]==0 & ... & CENS[i-1]==0)
    # (2) An alternative: collapse CENS into a categorical (automatically), based on the ordering in CENS. Then the fitting of categoricals will perform all the subsetting correctly.
    # (3) Alternative: set the indicators of missingness in the right way for CENS[i] if any CENS[1], ..., CENS[i-1] are 1.

  # - Stratification - allows K models on the SAME OUTCOME by stratifying rule
    # (1) User specified rule function creates strata. (stratify.CENS, stratify.TRT, stratify.MONITOR) Note that if the rule is based on data.table syntax it will be VERY FAST!
    # A nice trick would be to be able to AUTOMATICALLY convert logical subset expressions to data.table statements -> Its possible with some meta-programming and parsing
    # (2) These subsets (logical vectors) define K regressions, one regression model for each subset expression
    # Can specify K regressions in gform.CENS/gform.TRT/gform.MONITOR. If only one regression is specified it will be aplied to ALL stratas.

  # - Factors: when a predictor in covars is a factor it needs to be automatically converted to length(levels(L[i])) dummy indicators.
  # Might use the existing routine from tmlenet by factor -> categorical -> indicatorMatrix


# - consider new dcast.data.table features in v1.9.7. Could be useful for fast conversion to wide format -> use simcusal interface for creating summaries -> back to long format with melt.data.table.
 # now allows drop = c(FALSE, TRUE) and drop = c(TRUE, FALSE). The former only fills all missing combinations of formula LHS, where as the latter fills only all missing combinations of formula RHS.
# Thanks to Ananda Mahto for this SO post and to Jaap for filing #1512.
# http://stackoverflow.com/questions/34830908/make-the-drop-argument-in-dcast-only-look-at-the-rhs-of-the-formula
# https://github.com/Rdatatable/data.table
# https://github.com/Rdatatable/data.table/issues/1512

# - specifically look into g-force optimized functions for data.table: https://github.com/Rdatatable/data.table/issues/523


}












