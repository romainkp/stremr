#' Estimate the Survival of Intervention on Exposures and MONITORing Process for Right Censored Longitudinal Data.
#'
#' The \pkg{stremr} R package is a tool for estimation of causal survival curve under various user-specified interventions
#' (e.g., static, dynamic, deterministic, or stochastic).
#' In particular, the interventions may represent exposures to treatment regimens, the occurrence or non-occurrence of right-censoring
#' events, or of clinical monitoring events. \pkg{stremr} enables estimation of a selected set of the user-specified causal quantities of interest,
#' such as, treatment-specific survival curves and the average risk difference over time.
#'
#' @section Documentation:
#' \itemize{
#' \item To see the package vignette use: \code{vignette("stremr_vignette", package="stremr")}
#' \item To see all available package documentation use: \code{help(package = 'stremr')}
#' }
#'
#' @section Routines:
#' ...
#'
#' @section Data structures:
#' The following most common types of output are produced by the package:
#' \itemize{
#' \item \emph{observed data} - input data.frame in long format (repeated measures over time).
#' }
#'
#' @section Updates:
#' Check for updates and report bugs at \url{http://github.com/osofr/stremr}.
#'
#' @docType package
#' @name stremr-package
#'
NULL

#' An example of a dataset in long format with categorical censoring variable.
#'
#' Simulated dataset containing 1,000 i.i.d. observations organized in long format as person-time rows.
#' The binary exposure is \code{TI} and binary outcome is \code{Y.tplus1}. See /tests/ for R code that generated this data.
#'
#' @format A data frame with 1,000 observations and variables:
#' \describe{
#'   \item{ID}{Unique subject identifier}
#'   \item{CVD}{Baseline confounder (time invariant)}
#'   \item{t}{Interger for current time period, range 0-16}
#'   \item{lastNat1}{Time since last monitoring event, set to 0 when N[t-1]=0 and then added one for each new period where N[t] is 0.}
#'   \item{highA1c}{Time-varying confounder}
#'   \item{CatC}{Categorical censoring variable, range 0-2. The value of 0 indicates no censoring 1 or 2 indicates censoring (possibly for different reasons)}
#'   \item{C}{Binary censoring indicator derived from CatC. 0 if CatC is 0 and 1 if CatC is 1 or 2.}
#'   \item{TI}{Binary exposure variable}
#'   \item{N}{The indicator of being monitored (having a visit)}
#'   \item{Y.tplus1}{Binary outcome at t}
#' }
#' @docType data
#' @keywords datasets
#' @name OdataCatCENS
#' @usage data(OdataCatCENS)
NULL


#' An example of a dataset in long format with random monitoring and no right censoring.
#'
#' Simulated dataset containing 1,000 i.i.d. observations organized in long format as person-time rows.
#' The binary exposure is \code{TI} and binary outcome is \code{Y.tplus1}. See /tests/ for R code that generated this data.
#'
#' @format A data frame with 1,000 observations and variables:
#' \describe{
#'   \item{ID}{Unique subject identifier}
#'   \item{CVD}{Baseline confounder (time invariant)}
#'   \item{t}{Interger for current time period, range 0-16}
#'   \item{lastNat1}{Time since last monitoring event, set to 0 when N[t-1]=0 and then added one for each new period where N[t] is 0.}
#'   \item{highA1c}{Time-varying confounder}
#'   \item{TI}{Binary exposure variable}
#'   \item{C}{Administrative censoring indicator, always set to 0 unless the end of study is reached (t==16)}
#'   \item{N}{The indicator of being monitored (having a visit)}
#'   \item{Y.tplus1}{Binary outcome at t}
#' }
#' @docType data
#' @keywords datasets
#' @name OdataNoCENS
#' @usage data(OdataNoCENS)
NULL

#' An example of a dataset in long format with random monitoring process and no right censoring.
#'
#' Simulated dataset containing 10,000 i.i.d. observations organized in long format as person-time rows.
#' The binary exposure is \code{TI} and binary outcome is \code{Y.tplus1}. See /tests/
#' for R code that generated this data as well as R code that uses stremr to analyze this data.
#'
#' @format A data frame with 10,000 observations and variables:
#' \describe{
#'   \item{ID}{Unique subject identifier}
#'   \item{t}{Interger for current time period, range 0-16}
#'   \item{CVD}{Baseline confounder (time invariant)}
#'   \item{lastNat1}{Time since last monitoring event, set to 0 when N[t-1]=0 and then added one for each new period where N[t] is 0.}
#'   \item{highA1c}{Time-varying confounder}
#'   \item{TI}{Binary exposure variable}
#'   \item{C}{Administrative censoring indicator, always set to 0 unless the end of study is reached (t==16)}
#'   \item{N}{The random indicator of being monitored (having a visit), simulated as a Bernoulli RV with P(N(t)=1)=0.5}
#'   \item{Y.tplus1}{Indicator of the survival event at t}
#'   \item{gTI.dlow}{Counterfactual exposure under static intervention - always treat}
#'   \item{gTI.dhigh}{Counterfactual exposure under dynamic intervention - treat only when highA1c is above 1 and the subject is being monitored}
#'   \item{gPois3.yrly}{Poisson probability of counterfactual monitoring indicator being equal to 1}
#'   \item{gPois3.biyrly}{Poisson probability of counterfactual monitoring indicator being equal to 1}
#'   \item{gp05}{Bernoulli probability of counterfactual monitoring indicator being equal to 1}
#' }
#' @docType data
#' @keywords datasets
#' @name OdatDT_10K
#' @usage data(OdatDT_10K)
NULL

#' An example of a dataset in long format with two time-points and no censoring.
#'
#' Simulated dataset containing 500,000 i.i.d. observations organized in long format as person-time rows.
#' The binary time-varying exposure is \code{A} and binary time-invariant outcome is \code{Y}. See /tests/
#' for R code that generated this data as well as R code that uses stremr to analyze this data.
#'
#' @format A data frame with 500,000 observations and variables:
#' \describe{
#'   \item{ID}{Unique subject identifier}
#'   \item{t}{Interger for current time period, range 0-1}
#'   \item{L1}{Binary time-varying confounder}
#'   \item{L2}{Binary time-varying confounder}
#'   \item{L3}{Binary time-varying confounder}
#'   \item{A}{Binary time-varying exposure}
#'   \item{Y}{Indicator of the binary outcome at the end of the study (t=1). Constant set to 0 at first time point (t=0)}
#' }
#' @docType data
#' @keywords datasets
#' @name dt_2t
#' @usage data(dt_2t)
NULL


#' An example of a dataset in long format with one time-point and no censoring.
#'
#' Simulated dataset containing 50,000 i.i.d. observations with single time-point exposure.
#' The binary exposure is \code{A} and binary outcome is \code{Y}. See /tests/
#' for R code that generated this data as well as R code that uses stremr to analyze this data.
#'
#' @format A data frame with 50,000 observations and variables:
#' \describe{
#'   \item{ID}{Unique subject identifier}
#'   \item{L1}{Binary baseline confounder}
#'   \item{L2}{Binary baseline confounder}
#'   \item{L3}{Binary baseline confounder}
#'   \item{L4}{Binary baseline confounder}
#'   \item{L5}{Binary baseline confounder}
#'   \item{L6}{Binary baseline confounder}
#'   \item{A}{Binary exposure}
#'   \item{Y}{Indicator of the binary outcome at the end of the study.}
#' }
#' @docType data
#' @keywords datasets
#' @name dt_1t
#' @usage data(dt_1t)
NULL







