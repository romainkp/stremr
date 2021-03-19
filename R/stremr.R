#' @import R6
#' @import data.table
#' @importFrom graphics axis barplot hist par text  legend plot
#' @importFrom methods is
#' @importFrom stats approx binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str capture.output
#' @importFrom stats as.formula glm na.exclude rbinom terms.formula pnorm quasibinomial
NULL


if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("estimates", "FUPtimes_tabs", "MSM", "MSM.crude", "NPMSM", "trunc_MSM", "trunc_TMLE", "trunc_weight", "wts_data", "wts_tabs"))
}


tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
  invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler), warning = W)
}

#-----------------------------------------------------------------------------
# Class Membership Tests
#-----------------------------------------------------------------------------
is.DataStorageClass <- function(DataStorageClass) "DataStorageClass"%in%class(DataStorageClass)
is.ModelStack <- function(obj) {
  ("ModelStack" %in% class(obj)) || ("splitCVStack" %in% class(obj))
}

#-----------------------------------------------------------------------------
# Capture the arguments passed on as ... in a list
#-----------------------------------------------------------------------------
capture.exprs <- function(...) {
  # sVar.exprs <- eval(substitute(alist(...)))
  sVar.exprs <- list(...)
  if (is.null(names(sVar.exprs))) names(sVar.exprs) <- rep_len("", length(sVar.exprs))
  if (length(sVar.exprs)!=0 && any(names(sVar.exprs)%in%"")) {
    stop("all parameters passed to ... must be named")
  }
  return(sVar.exprs)
}

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
# MAIN stremr FUNCTION
#---------------------------------------------------------------------------------
#' Run all estimators (KM, hazard-IPW, direct-IPW, GCOMP plug-in and TMLE).
#'
#' Estimate the causal survival curve for a particular stochastic, dynamic or static intervention on the treatment/exposure and monitoring process.
#'  Implements the \strong{IPW} (Inverse Probability-Weighted or Horvitz-Thompson) estimator of the discrete survival hazard function which is mapped into survival function.
#' @param data Input data in long format. Can be a \code{data.frame} or a \code{data.table} with named columns, containing the time-varying covariates (\code{covars}),
#'  the right-censoring event indicator(s) (\code{CENS}), the exposure variable(s) (\code{TRT}), the monitoring process variable(s) (\code{MONITOR})
#'  and the survival OUTCOME variable (\code{OUTCOME}).
# @param estimators (NOT IMPLEMENTED) Character vector with estimator names.
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
#' @param OUTCOME  A column name in \code{data} for the survival OUTCOME variable name, code as 1 for the outcome event.
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
# @param intervened_MONITOR Column name in \code{data} containing the counterfactual probabilities of following a specific monitoring regimen.
#' @param noCENScat Same as in \code{\link{importData}}.
#' @param remove_extra_rows Same as in \code{\link{importData}}.
#' @param nfolds Number of folds for cross-validation (leave as \code{NULL} if no cross-validation is desired).
#' @param models_CENS Same as in \code{\link{fitPropensity}}.
#' @param models_TRT Same as in \code{\link{fitPropensity}}.
#' @param models_MONITOR Same as in \code{\link{fitPropensity}}.
#' @param fit_method_g Same as \code{fit_method} in \code{\link{fitPropensity}}.
#' @param models_Q Same as \code{models} in \code{\link{fit_GCOMP}}.
#' @param fit_method_Q Same as \code{fit_method} in \code{\link{fit_GCOMP}}.
#' @param Qforms Same as in \code{\link{fit_GCOMP}}.
#' @param tvals Same as in \code{\link{fit_GCOMP}}.
#' @param stratifyQ_by_rule Same as in \code{\link{fit_GCOMP}}.
#' @param trunc_IPW_MSM Weight truncation for IPW-based functions.
#' @param trunc_IPW_TMLE Weight trunction for TMLE.
#' @param seed Random generator seed.
#' @param MSMGLMpkg Package to use for MSM GLM fits (see \code{\link{survMSM}}).
#' @param tbreaks Same as in \code{\link{survMSM}}.
#' @param start_h2o_cluster Start h2o cluster?
#' @param nthreads Number of threads (CPUs) to use for h2o cluster?
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(stremr.verbose=TRUE)}.
#' @return ...
#' @seealso \code{\link{stremr-package}} for the general overview of the package,
#' @example tests/examples/1_stremr_example.R
#' @export
stremr <- function(data, ID = "Subj_ID", t_name = "time_period",
                   covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
                   gform_CENS, gform_TRT, gform_MONITOR,
                   stratify_CENS = NULL, stratify_TRT = NULL, stratify_MONITOR = NULL,
                   intervened_TRT = NULL,
                   # intervened_MONITOR = NULL,
                   noCENScat = 0L,
                   remove_extra_rows = TRUE,
                   nfolds = NULL,
                   models_CENS = sl3::Lrnr_glm_fast$new(family = quasibinomial()),
                   # models_CENS = gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial"),
                   models_TRT = sl3::Lrnr_glm_fast$new(family = quasibinomial()),
                   # models_TRT = gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial"),
                   models_MONITOR = sl3::Lrnr_glm_fast$new(family = quasibinomial()),
                   # models_MONITOR = gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial"),
                   fit_method_g = "none",
                   models_Q = sl3::Lrnr_glm_fast$new(family = quasibinomial()),
                   # models_Q = gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial"),
                   fit_method_Q = "none",
                   Qforms = NULL,
                   tvals = NULL,
                   stratifyQ_by_rule = TRUE,
                   trunc_IPW_MSM = Inf,
                   trunc_IPW_TMLE = Inf,
                   seed = NULL,
                   MSMGLMpkg = c("speedglm", "h2o"),
                   tbreaks = NULL,
                   start_h2o_cluster = TRUE,
                   nthreads = 2,
                   verbose = getOption("stremr.verbose")) {

  if (start_h2o_cluster)
    h2o::h2o.init(nthreads = nthreads)


  ## regression formulas for Q's:
  # Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(tvals)+1))
  if (is.null(Qforms)) {
    Qforms <- "Qkplus1 ~ " %+% paste0(c(covars, CENS, TRT, MONITOR), collapse = " + ")
    message("Using default Qforms: " %+% Qforms)
  }

  # ------------------------------------------------------------------
  # - BLOCK 1: Process inputs and define OData R6 object
  # ------------------------------------------------------------------
  OData <- importData(data, ID = ID, t_name = t_name, covars = covars,
                      CENS = CENS, TRT = TRT, MONITOR = MONITOR,
                      OUTCOME = OUTCOME, noCENScat = noCENScat,
                      remove_extra_rows = remove_extra_rows, verbose = verbose)

  if (!is.null(nfolds)) {
    assert_that(is.integer(nfolds))
    OData <- define_CVfolds(OData, nfolds = nfolds, fold_column = "fold_ID", seed = seed)
  }

  if (is.null(tvals)) {
    tvals <- OData$min.t:OData$max.t
    message("Using default tvals: " %+% paste(tvals, collapse=","))
  }

  if (is.null(tbreaks)) {
      tbreaks <- OData$min.t:(OData$max.t-1)
      message("Using default tbreaks: " %+% paste(tbreaks, collapse=","))
    }

  # ------------------------------------------------------------------
  # - BLOCK 2: define regression models, define a single RegressionClass & fit the propensity score for observed data, summary.g0 g0 (C,A,N)
  # ------------------------------------------------------------------
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR,
                         stratify_CENS = stratify_CENS, stratify_TRT = stratify_TRT, stratify_MONITOR = stratify_MONITOR,
                         models_CENS = models_CENS, models_TRT = models_TRT, models_MONITOR = models_MONITOR,
                         fit_method = fit_method_g)
                          # ,
                          # fold_column = fold_column)

  ## ---------------------------------------------------------------------------------------
  ## - BLOCK 3: Define analyses data set
  ## **** As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** This dataset is to be saved and will be later merged in with all analysis
  ## ---------------------------------------------------------------------------------------
  ## This dataset defines all parameters that we like to vary in this analysis (including different interventions)
  ## That is, each row of this dataset corresponds with a single analysis, for one intervention of interest.

  analysis <- list(intervened_TRT = intervened_TRT,
                   # intervened_MONITOR = intervened_MONITOR,
                  stratifyQ_by_rule = stratifyQ_by_rule) %>%
                  purrr::cross_df() %>%
                  dplyr::arrange(stratifyQ_by_rule) %>%
                  dplyr::mutate(trunc_MSM = trunc_IPW_MSM[1L]) %>%
                  dplyr::mutate(trunc_TMLE = trunc_IPW_TMLE[1L])

  ## ---------------------------------------------------------------------------------------
  ## - BLOCK 4:
  ## Evaluate weights based gstar_TRT, gstar_MONITOR and observed propensity scores g0, the input is modelfits.g0 and OData object
  ## Non-parametric MSM for survival, with weight stabilization, input either single weights dataset or a list of weights datasets,
  ## Each dataset containing weights non-zero weights for single regimen
  ## ---------------------------------------------------------------------------------------
  IPW <- analysis %>%
    dplyr::rename(trunc_weight = trunc_MSM) %>%
    dplyr::distinct(intervened_TRT, trunc_weight) %>%
    dplyr::group_by(intervened_TRT) %>%
    dplyr::mutate(wts_data = purrr::map(dplyr::first(intervened_TRT), getIPWeights, OData = OData, tmax = OData$max.t)) %>%
    ## save the tables of weights summaries (sep for each regimen)
    dplyr::mutate(wts_tabs = purrr::map(wts_data,
        ~ get_wtsummary(.x, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE))) %>%
    ## save the tables with number at risk / following each rule (sep for each regimen)
    dplyr::mutate(FUPtimes_tabs = purrr::map(wts_data,
          ~ get_FUPtimes(.x, IDnode = ID, tnode = t_name))) %>%
    dplyr::ungroup() %>%

    ## IPW-Adjusted KM (Non-Parametric or Saturated MSM):
    dplyr::mutate(NPMSM = purrr::map2(wts_data, trunc_weight,
      ~ survNPMSM(wts_data = .x,
                  trunc_weights = .y,
                  OData = OData))) %>%
    dplyr::mutate(NPMSM = purrr::map(NPMSM, "estimates")) %>%

    ## IPW-Adjusted KM (Non-Parametric or Saturated MSM):
    dplyr::mutate(NPMSM = purrr::map2(wts_data, trunc_weight,
      ~ survNPMSM(wts_data = .x,
                  trunc_weights = .y,
                  OData = OData))) %>%
    dplyr::mutate(NPMSM = purrr::map(NPMSM, "estimates")) %>%

    ## Crude MSM for hazard (w/out IPW):
    dplyr::mutate(MSM.crude = purrr::map(wts_data,
      ~ survMSM(wts_data = .x,
                OData = OData,
                tbreaks = tbreaks,
                use_weights = FALSE,
                glm_package = MSMGLMpkg[1L]))) %>%
    dplyr::mutate(MSM.crude = purrr::map(MSM.crude, "estimates")) %>%

    ## IPW-MSM for hazard (smoothing over time-intervals in tbreaks):
    dplyr::mutate(MSM = purrr::map2(wts_data, trunc_weight,
      ~ survMSM(wts_data = .x,
                trunc_weights = .y,
                OData = OData,
                tbreaks = tbreaks,
                glm_package = MSMGLMpkg[1L]))) %>%
    dplyr::mutate(MSM = purrr::map(MSM, "estimates")) %>%

    dplyr::mutate(directIPW = purrr::map2(wts_data, trunc_weight,
      ~ directIPW(wts_data = .x,
                  trunc_weights = .y,
                  OData = OData))) %>%
    dplyr::mutate(directIPW = purrr::map(directIPW, "estimates")) %>%
    dplyr::rename(trunc_MSM = trunc_weight)

    ## save IPW tables (will be later merged with main results dataset)
    IPWtabs <-  analysis %>%
      dplyr::left_join(IPW) %>%
      dplyr::distinct(intervened_TRT, trunc_MSM, wts_tabs, FUPtimes_tabs) %>%
      tidyr::nest(intervened_TRT, wts_tabs, FUPtimes_tabs, .key = "IPWtabs")

    IPW <- IPW %>% dplyr::select(-wts_data, -wts_tabs, -FUPtimes_tabs)

  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  if (length(Qforms)==1) {
    Qforms <- rep.int(Qforms, length(tvals))
  }

  GCOMP <-analysis %>%
    dplyr::distinct(intervened_TRT, stratifyQ_by_rule) %>%
    dplyr::mutate(GCOMP = purrr::map2(intervened_TRT, stratifyQ_by_rule,
          ~ fit_GCOMP(intervened_TRT = .x,
                     stratifyQ_by_rule = .y,
                     tvals = tvals,
                     OData = OData,
                     models = models_Q,
                     Qforms = Qforms,
                     fit_method = fit_method_Q))) %>%
    dplyr::mutate(GCOMP = purrr::map(GCOMP, "estimates"))

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  TMLE <- analysis %>%
    dplyr::rename(trunc_weight = trunc_TMLE) %>%
    dplyr::distinct(intervened_TRT, stratifyQ_by_rule, trunc_weight)

  TMLE <- TMLE %>%
    dplyr::mutate(TMLE = purrr::pmap(TMLE, fit_TMLE,
                   tvals = tvals,
                   OData = OData,
                   models = models_Q,
                   Qforms = Qforms,
                   fit_method = fit_method_Q)) %>%
    dplyr::mutate(TMLE = purrr::map(TMLE, "estimates")) %>%
    dplyr::rename(trunc_TMLE = trunc_weight)

  ## ------------------------------------------------------------
  ## COMBINE ALL ANALYSES INTO A SINGLE DATASET
  ## ------------------------------------------------------------
  results <-  analysis %>%
              dplyr::left_join(IPW) %>%
              dplyr::left_join(GCOMP) %>%
              dplyr::left_join(TMLE)

  ## Nest each estimator by treatment regimen (we now only show the main analysis rows)
  results <- results %>%
             tidyr::nest(intervened_TRT, NPMSM, MSM.crude, MSM, directIPW, GCOMP, TMLE, .key = "estimates")

  ## Calculate RDs (contrasting all interventions, for each analysis row & estimator).
  ## The RDs data no longer needs the intervened_TRT column
  results <-  results %>%
              dplyr::mutate(RDs =
                purrr::map(estimates,
                  ~ dplyr::select(.x, -intervened_TRT) %>%
                  purrr::map(~ get_RDs(.x)) %>%
                  tibble::as_tibble()
                  ))

  # ---------------------------------------------------------------------------------------
  # - BLOCK 5: Builds a report with weight distributions, survival estimates, etc.
  # ---------------------------------------------------------------------------------------

return(results)
}
