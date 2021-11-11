#' @importFrom magrittr %>%
NULL

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "dx1", "dx2", "pval", "CI", "estimator"))
}

#' Evaluate the long format dataset of risk differences over estimated survival time-points
#'
#' The risk difference (RD) of regimens 1 (dx1) and 2 (dx2) at time-point \code{t} is
#' calculated as \code{S_2(t) - S_1(t)},
#' where \code{S_2(t)} and \code{S_1(t)} are the corresponding survival estimates
#' at \code{t} for regimens 2 and 1, respectively.
#' The RDs are evaluated for all intersecting time-points in the input data
#' \code{St_data}. When the input data contains survival estimates for more
#' than 2 regimens the RDs between all regimens are evaluated in their order
#' of appearance. This order can be controlled with the argument \code{order}.
#' @param St_data A list containing \code{data.table} survival estimates for
#' each regimen to be contrasted.
#' @param St_name The name of the column containing the survival estimates.
#' @param getSEs Should standard errors of risk difference also be estimates?
#' Note that the estimates will only be available when the data.table for
#' each regimen contains a column of subject-specific influence curve (IC)
#' based estimates. The column containing these ICs should be named \code{"IC.St"}.
#' @param order In which order should the RD contrasts be evaluated?
#' For example, suppose the input data contains
#' the survival estimates for 3 regimens, where \code{dx1}, \code{dx2}, \code{dx3}
#' denote these three regimens and \code{S_1(t)}, \code{S_2(t)}, \code{S_3(t)}
#' denote the corresponding  survival estimates at time-point \code{t}.
#' By default, the following three contrasts will be evaluated in this order:
#' RD(dx1,dx2) = \code{S_2(t)-S_1(t)}, RD(dx1,dx3)=\code{S_3(t)-S_1(t)} and
#' RD(dx2,dx3)=\code{S_3(t)-S_2(t)}.
#' However, the argument \code{order} can be used to change which of the possible
#' three contrasts are evaluated.
#' For example, when setting \code{order=c(3,2,1)}, the following 3 RDs are evaluated instead:
#' RD(dx3,dx2), RD(dx3,dx1) and RD(dx2,dx1).
#' @return A long format dataset (\code{tibble}) with risk differences.
#' @export
get_RDs <- function(St_data, St_name, getSEs = TRUE, order = seq_along(St_data)) {

  ## Use default name for the column with survival estimates
  if (missing(St_name)) {
    estimator_short <- unique(attr(St_data[[1]], "estimator_short"))
    if (is.null(estimator_short)) estimator_short <- unique(St_data[[1]][["est_name"]])
    if (is.null(estimator_short)) {
      stop("Cannot automatically detect the name of the column that contains the estimates. Please provide the name using 'St_name = ...'.")
    }
    St_name <- paste0("St.", estimator_short)
  }

  ## check that the influence curve column exists in all survival datasets
  ## if not, automatically set getSEs to FALSE
  if (!all(unlist(purrr::map( St_data, ~ "IC.St" %in% names(.x) ))))
    getSEs <- FALSE

  # nIDs <- attr(St_data[[1]], "nID")

  ## Grab all available follow-up time-points by regimen:
  # time_bydx_obs <- lapply(St_data, function(St) attr(St, "time"))
  time_bydx_obs <- lapply(St_data, function(St) St[["time"]])

  ## Find the intersection of follow-up time-points across all regimens
  ## Will evaluate all RDs only at these intersected f-up time-points
  time <- Reduce(intersect, time_bydx_obs)
  time_idx <- seq_along(time)

  if (!all(unlist(lapply(St_data, function(surv) St_name %in% names(surv)))))
    stop("name of the survival estimates column cannot be found in one of the input datasets: " %+% St_name)

  tx_idx <- seq_along(St_data)
  # tx_names <- unlist(lapply(St_data, function(St) attr(St, "rule_name")))
  tx_names <- unlist(lapply(St_data, function(St) unique(St[["rule.name"]])))

  ## TO EVALUATE RD (dx1 - dx2) IS THE SAME AS EVALUATING S.t_dx2 - S.t_dx1
  eval_RDs_two_tx <- function(dx1, dx2, time_idx, St_name, ...) {
    St_data[[dx2]][[St_name]][time_idx] - St_data[[dx1]][[St_name]][time_idx]
  }

  eval_SEs_two_tx <- function(dx1, dx2, time_idx, ...) {
    nIDs_1 <- length(St_data[[dx1]][["IC.St"]][[time_idx]])
    nIDs_2 <- length(St_data[[dx2]][["IC.St"]][[time_idx]])
    if (nIDs_1 != nIDs_2) {
      warning("Cannot evaluate SEs for RDs from two ICs for dx1='" %+% dx1 %+% "' and dx2='" %+% dx2 %+%"', since their lengths differ: " %+% nIDs_1 %+% " vs. " %+% nIDs_2 %+% ". 
  Evaluating only the point RD estimates.")
      return(NA_real_)
    } else {
      nIDs <- nIDs_1
      SE.diff <- sqrt(sum(( St_data[[dx2]][["IC.St"]][[time_idx]] - St_data[[dx1]][["IC.St"]][[time_idx]] )^2) / nIDs^2)
      return(SE.diff)
    }
  }

  gs <- list(dx1 = tx_idx,
             dx2 = tx_idx,
             time_idx = time_idx) %>%
        purrr::cross_df() %>%
        dplyr::mutate(rank1 = order[dx1]) %>%
        dplyr::mutate(rank2 = order[dx2]) %>%
        dplyr::filter(rank1 < rank2) %>%
        dplyr::arrange(rank1, rank2) %>%
        dplyr::select(-rank1, -rank2) %>%
        dplyr::mutate(time = time[time_idx])

  gs <- gs %>%
        dplyr::mutate(RD = purrr::pmap_dbl(., eval_RDs_two_tx, St_name)) %>%
        dplyr::mutate(RD.SE = NA_real_)

  if (getSEs)
    gs <- gs %>%
          dplyr::mutate(RD.SE = purrr::pmap_dbl(., eval_SEs_two_tx))

    gs <- gs %>%
          dplyr::mutate(CI95low = round(RD - abs(qnorm(0.025))*RD.SE, 4)) %>%
          dplyr::mutate(CI95up = round(RD + abs(qnorm(0.025))*RD.SE, 4)) %>%
          tidyr::unite("CI", CI95low, CI95up, sep = ";", remove = FALSE) %>%
          dplyr::mutate(pval = round(2*pnorm( abs(RD/RD.SE), lower.tail=FALSE ),2))

  gs <- gs %>%
        dplyr::mutate(dx1_name = purrr::map_chr(dx1, ~ tx_names[.x])) %>%
        dplyr::mutate(dx2_name = purrr::map_chr(dx2, ~ tx_names[.x])) %>%
        tidyr::unite("contrast", dx1_name, dx2_name, remove = FALSE) %>%
        dplyr::mutate(estimator = St_name)

  attr(gs, "tx") <- tx_names
  attr(gs, "stratifyQ_by_rule") <- attr(St_data[[1]], "stratifyQ_by_rule")
  attr(gs, "trunc_weights") <- attr(St_data[[1]], "trunc_weights")

  return(list(gs))
}

#' Print RD table
#'
#' Print and return RD table in various formats using \code{knitr::kable}
#' @param RD_tab Risk difference table returned by \code{get_RDs}.
#' @param dx1 Subset \code{RD_tab} by these \code{dx1} values.
#' @param dx2 Subset \code{RD_tab} by these \code{dx2} values.
#' @param time Subset \code{RD_tab} by these \code{time} values.
#' @param ... Additional parameters to be passed directly to \code{knitr::kable}.
#' @return A table returned by calling \code{knitr::kable()}.
#' @export
print_RDs <- function(RD_tab, dx1, dx2, time, ...) {
  if (missing(dx1)) {
    dx1_sel <- unique(RD_tab[["dx1"]])
  } else {
    dx1_sel <- dx1
  }

  if (missing(dx2)) {
    dx2_sel <- unique(RD_tab[["dx2"]])
  } else {
    dx2_sel <- dx2
  }

  if (missing(time)) {
    time_sel <- unique(RD_tab[["time"]])
  } else {
    time_sel <- time
  }

  RD_tab_sel <- 
    RD_tab %>% 
    dplyr::filter((dx1 %in% dx1_sel) & (dx2 %in% dx2_sel) & (time %in% time_sel)) %>%
    dplyr::select(estimator, dx1_name, dx2_name, time, RD, RD.SE, CI, pval)

  var <- knitr::kable(RD_tab_sel, ...)
  print(var)
  return(var)
}


#' Follow-up times by regimen
#'
#' Subject specific follow-up times for each regimen in \code{wts_data}.
#' @param wts_data Either a list of data.table containing weights (one for each separate regimen/intervention) or a single data.table with
#' weights for one regimen / intervention.
#' @param IDnode Name of the column containing subject-specific identifier in the input data.
#' @param tnode Name of the column containing the time/period variable in the input data.
#' @return A \code{data.table} of subject specific follow-up times for each regimen in \code{wts_data}.
#' @seealso \code{\link{getIPWeights}} for evaluation of IP-weights.
#' @export
get_FUPtimes <- function(wts_data, IDnode, tnode) {
  wts_data <- format_wts_data(wts_data)
  t.name.col <- tnode
  ID.name.col <- IDnode
  follow_up_rule_ID <- wts_data[cum.IPAW > 0, list(max.t = max(get(t.name.col), na.rm = TRUE)), by = list(get(ID.name.col), get("rule.name"))]
  data.table::setnames(follow_up_rule_ID, c(IDnode, "rule.name", "max.t"))
  data.table::setkeyv(follow_up_rule_ID, cols = IDnode)
  return(follow_up_rule_ID)
}

#' IP-Weights Summary Tables
#'
#' Produces various table summaries of IP-Weights.
#' @param wts_data Either a list of data.table containing weights (one for each separate regimen/intervention) or a single data.table with
#' weights for one regimen / intervention.
#' @param cutoffs Weight cut off points for summary tables.
#' @param varname Character string describing the type of the weights
#' @param by.rule Can optionally evaluate the same summary tables separately for each regimen / rule.
#' @param stabilize Set to \code{TRUE} to return stabilized weights summary, otherwise unstabilized weights are used (default).
#' @return A list with various IP-weights summary tables.
#' @seealso \code{\link{getIPWeights}} for evaluation of IP-weights.
#' @export
get_wtsummary <- function(wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), varname = "Stabilized IPAW", by.rule = FALSE, stabilize = FALSE) {
  wts_data_new <- copy(wts_data)
  wts_data_new <- format_wts_data(wts_data_new)
  if (stabilize) wts_data_new[, `:=`("cum.IPAW", cum.stab.P * cum.IPAW)]
  # get the counts for each category (bin) + missing count:
  na.count <- sum(is.na(wts_data_new[["cum.IPAW"]]))
  na.yes <- (na.count > 0L)
  # define a list of intervals:
  cutoffs2 <- c(-Inf, cutoffs, Inf)
  idx <- seq_along(cutoffs2); idx <- idx[-length(idx)]
  intervals_list <- lapply(idx, function(i) c(cutoffs2[i], cutoffs2[i+1]))
  x.freq.counts.DT <- wts_data_new[, lapply(intervals_list, function(int) sum((cum.IPAW >= int[1]) & (cum.IPAW < int[2]), na.rm = TRUE))]
  x.freq.counts <- as.integer(x.freq.counts.DT[1,])
  if (na.yes) {
    x.freq.counts <- c(x.freq.counts, na.count)
    x.freq.counts.DT <- cbind(x.freq.counts.DT, wts_data_new[, sum(is.na(cum.IPAW), na.rm = TRUE)])
  }
  # create labels:
  catNames <- rep.int(NA, (length(cutoffs) + 1 + as.integer(na.yes)))
  catNames[1] <- paste("<", cutoffs[1], sep="")
  for(i in 1:(length(cutoffs)-1)){
    catNames[i+1] <- paste("$[$",cutoffs[i],", ",cutoffs[i+1],"$[$",sep="")
  }
  catNames[length(cutoffs)+1] <- paste(">=",cutoffs[length(cutoffs)],sep="")
  if (na.yes) catNames[length(cutoffs) + 2] <- "Missing"
  names(x.freq.counts) <- catNames
  colnames(x.freq.counts.DT) <- catNames
  # make a table:
  x.freq.counts.tab <- makeFreqTable(x.freq.counts)
  colnames(x.freq.counts.tab)[1] <- varname
  x.freq.counts.tab[length(cutoffs)+1,1] <- paste("$\\geq$ ",cutoffs[length(cutoffs)],sep="")
  # do the same by separately for each rule:
  if (by.rule) {
    x.freq.counts.byrule.DT <- wts_data_new[, lapply(intervals_list, function(int) sum((cum.IPAW >= int[1]) & (cum.IPAW < int[2]), na.rm = TRUE)), by = rule.name]
    if (na.yes) {
      x.freq.counts.byrule.DT <- x.freq.counts.byrule.DT[wts_data_new[, sum(is.na(cum.IPAW), na.rm = TRUE), by = rule.name], on = c("rule.name")]
    }
    colnames(x.freq.counts.byrule.DT)[-1] <- catNames
  } else {
    x.freq.counts.byrule.DT <- NULL
  }
  return(list(summary.table = x.freq.counts.tab, summary.DT = x.freq.counts.DT, summary.DT.byrule = x.freq.counts.byrule.DT))
}











# ---------------------------------------------------------------------------------------------
#' Helper routine to define the monitoring indicator and time since last visit
#'
#' @param data Input \code{data.table} or \code{data.frame}
#' @param ID The name of the unique subject identifier (character, numeric or factor).
#' @param t The name of the time/period variable in \code{data}.
#' @param I The name of the numeric biomarker value which determines the dynamic treatment rule at each time point t.
#' @param imp.I The name of the binary indicator of missingness or imputation for I at time point t. This is used for coding MONITOR(t-1):=1-imp.I(t).
#'  When imp.I(t)=1 it means that the patient was not observed (no office visit) at time-point t and hence no biomarker was measured.
#' @param MONITOR.name The name of the new column that represents for each row t the indicator of being MONITORed (having a visit) at time points t+1.
#' This variable is added as a new column to the output dataset.
#' The column MONITOR(t) is set to 1 when the indicator imp.I (imputed biomarker) is 0 at time-point t+1 and vice versa.
#' @param tsinceNis1 The name of the future column (created by this routine) that counts number of periods since last monitoring event at t-1. More precisely,
#' it is a function of the past \code{N(t-1)}, where 0 means that N(t-1)=1; 1 means that N(t-2)=1 and N(t-1)=0; etc.
#'
#' @section Details:
#'
#' Convert the input long format data with the time ordering: (I(t), imp.I(t), C(t), A(t))
#' into the data format required by the stremr input functions: (I(t), C(t), A(t), N(t):=1-imp.I(t+1)).
#' N(t) at time-point t is defined as the indicator of being observed (having an office visit) at time point t+1 (next timepoint after t)
#' The very first value of I(t) (at the first time-cycle) is ALWAYS ASSUMED observed/measured (hence I.imp=0 for each first subject-time observation).
#'
#' @section Data Format:
#'
#' The input data.frame data needs to be in long format.
#'
#' The format of the specified columns needs to be as follows.
#'
#' The time ordering of the input data at each t is as follows: (I, imp.I, CENS, TRT)
#'
#' The time ordering of the output data at each t is as follows: (I, CENS, TRT, MONITOR), where MONITOR(t)=1-imp.I(t+1).
#'
#' In output data.table, MONITOR(t-1)=1 indicates that the biomarker I(t) at t is observed and vice versa.
#'
#' In addition the output data.table will contain a column 'tsinceNis1', where:
#' \itemize{
#'   \item tsinceNis1(t) = 0 means that the person was monitored at time-point t-1.
#'   \item tsinceNis1(t) > 0 is the count of the number of cycles since last monitoring event.
#' }
#' @return A \code{data.table} in long format with ordering (I, CENS, TRT, MONITOR)
#' @export
defineMONITORvars <- function(data, ID, t, imp.I, MONITOR.name = 'N', tsinceNis1 = 'last.Nt'){
  ID.expression <- as.name(ID)
  indx <- as.name("indx")
  if (is.data.table(data)) {
    DT <- data.table(data, key=c(ID, t))
  } else if (is.data.frame(data)) {
    DT <- data.table(data, key=c(ID, t))
  } else {
    stop("input data must be either a data.table or a data.frame")
  }
  CheckVarNameExists(DT, ID)
  CheckVarNameExists(DT, t)
  CheckVarNameExists(DT, imp.I)
  # "Leading" (shifting up) and inverting indicator of observing I, renaming it to MONITOR value;
  # N(t-1)=1 indicates that I(t) is observed. Note that the very first I(t) is assumed to be always observed.
  DT[, (MONITOR.name) := shift(.SD, n=1L, fill=NA, type="lead"), by = eval(ID.expression), .SDcols = (imp.I)]
  DT[, (MONITOR.name) := 1L - get(MONITOR.name)]

  # Create "indx" vector that goes up by 1 every time MONITOR.name(t-1) shifts from 1 to 0 or from 0 to 1
  DT[, ("indx") := cumsum(c(FALSE, (get(MONITOR.name)!=0L)[-.N])), by = eval(ID.expression)]
  ## This line had a bug, since vec=cumsum(0, MONITOR) evaluates to length [.N+1].
  ## By applying vec[-.N] we remove the second to last observation, rather than just the last observation in vec
  # DT[, ("indx") := cumsum(c(FALSE, get(MONITOR.name)!=0L))[-.N], by = eval(ID.expression)]

  ## Evaluate the number of time-points that have passed since last visit:
  DT[, (tsinceNis1) := seq(.N)-1, by = list(eval(ID.expression), eval(indx))]
  DT[is.na(DT[["indx"]]), (tsinceNis1) := NA]
  DT[, ("indx") := NULL]
  return(DT)
}

# ---------------------------------------------------------------------------------------------
#' Define counterfactual dynamic treatment
#'
#' Defines a new column that contains the counterfactual dynamic treatment values.
#' Subject is switched to treatment when the biomarker \code{I} crosses the threshold \code{theta} while being monitored \code{MONITOR(t-1)=1}.
#' @param data Input data.frame or data.table in long format, see below for the description of the assumed format.
#' @param theta The vector of continuous cutoff values that index each dynamic treatment rule
#' @param ID The name of the unique subject identifier
#' @param t The name of the variable indicating time-period
#' @param I Continuous biomarker variable used for determining the treatmet decision rule
#' @param CENS Binary indicator of being censored at t;
#' @param TRT Binary indicator of the treatment (exposure) at t;
#' @param MONITOR The indicator of having a visit and having measured or observed biomarker I(t+1) (the biomarker value at THE NEXT TIME CYCLE).
#'  In other words the value of MONITOR(t-1) (at t-1) being 1 indicates that I(t) at time point t was observed/measured.
#'  The very first value of I(t) (at the first time-cycle) is ALWAYS ASSUMED observed/measured.
#' @param tsinceNis1 Character vector for the column in data, same meaning as described in convertdata(), must be already defined
#' @param new.TRT.names Vector with names which will be assigned to new columns generated by this routine(must be the same dimension as theta).
#'  When not supplied the following convention is adopted for naming these columns: \code{paste0(TRT,".gstar.","d",theta)}.
#' @param return.allcolumns Set to \code{TRUE} to return the original data columns along with new columns that define each rule
#' (can be useful when employing piping/sequencing operators).
#'
#' @section Details:
#'
#' * This function takes an input \code{data.frame} or \code{data.table}
#'   and produces an output data.table with counterfactual treatment assignment based on a rule defined by values of input column \code{I} and an input scalar \code{theta}.
#'
#' * Evaluates which observations should have received (switched to) treatment based on dynamic-decision rule defined by the measured biomarker \code{I}
#'   and pre-defined cutoffs supplied in the vector \code{theta}.
#'
#' * Produces a separate column for each value in \code{theta}:
#'\itemize{
#' \item (1) By default no one is treated (TRT assigned to 0)
#' \item (2) Subject may switch to treatment the first time \code{I} goes over the threshold \code{theta} (while having been monitored, i.e, MONITOR(t-1)==1)
#' \item (3) Once the subject switches to treatment he or she has to continue receiving treatment until the end of the follow-up
#'}
#'
#' * The format (time-ordering) of data is the same as required by the stremr() function: (I(t), CENS(t), TRT(t), MONITOR(t)).
#'   MONITOR(t) at time-point t is defined as the indicator of being observed (having an office visit) at time point t+1 (next timepoint after t)
#'   It is assumed that MONITOR(t) is always 0 for the very first time-point.
#'
#' @examples
#'
#' \dontrun{
#' theta <- seq(7,8.5,by=0.5)
#' FOLLOW.D.DT <- defineIntervedTRT(data = data, theta = theta,
#'                  ID = "StudyID", t = "X_intnum", I = "X_a1c",
#'                  TRT = "X_exposure", CENS = "X_censor", MONITOR = "N.t",
#'                  new.TRT.names = paste0("gstar",theta))
#' }
#' @example tests/examples/0_defineIntervention.R
#' @return A data.table with a separate column for each value in \code{theta}. Each column consists of indicators of following/not-following
#'  each rule indexed by a value form \code{theta}. In addition, the returned data.table contains \code{ID} and \code{t} columns for easy merging
#'  with the original data.
#' @export
defineIntervedTRT <- function(data, theta, ID, t, I, CENS, TRT, MONITOR, tsinceNis1, new.TRT.names = NULL, return.allcolumns = FALSE){
  ID.expression <- as.name(ID)
  chgTRT <- as.name("chgTRT")
  if (return.allcolumns) {
    DT <- data.table(data, key=c(ID,t))
  } else if (is.data.table(data)) {
    DT <- data.table(data[,c(ID, t, TRT, CENS, I, MONITOR, tsinceNis1), with = FALSE], key=c(ID,t))
  } else if (is.data.frame(data)) {
    DT <- data.table(data[,c(ID, t, TRT, CENS, I, MONITOR, tsinceNis1)], key=c(ID,t))
  } else {
    stop("input data must be either a data.table or a data.frame")
  }
  CheckVarNameExists(DT, ID)
  CheckVarNameExists(DT, t)
  CheckVarNameExists(DT, I)
  CheckVarNameExists(DT, TRT)
  CheckVarNameExists(DT, MONITOR)
  CheckVarNameExists(DT, tsinceNis1)
  if (!is.null(new.TRT.names)) {
    stopifnot(length(new.TRT.names)==length(theta))
  } else {
    new.TRT.names <- paste0(TRT,".gstar.","d",theta)
  }
  for (dtheta.i in seq_along(theta)) {
    dtheta <- theta[dtheta.i]
    # 1. By default, nobody is treated
    DT[, "new.TRT.gstar" :=  NA_integer_]
    # 2. The person can switch to treatment the earliest time 'I' goes over the threshold (while being monitored)
    # NOTE: ****** (tsinceNis1 > 0) is equivalent to (N(t-1)==0) ******
    DT[(get(tsinceNis1) == 0L) & (get(I) >= eval(dtheta)), "new.TRT.gstar" :=  1L]
    # 3. Once the person goes on treatment he/she has to stay on it until the end of the follow-up (using carry-forward)

    DT[, ("new.TRT.gstar") := zoo::na.locf(new.TRT.gstar, na.rm = FALSE), by = eval(ID.expression)]
    # DT[, new.TRT.gstar := zoo::na.locf(new.TRT.gstar, na.rm = FALSE), by = eval(ID.expression)]

    # 4. all remaining NA's must be the ones that occurred prior to treatment switch -> all must be 0 (not-treated)
    DT[is.na(new.TRT.gstar), ("new.TRT.gstar") := 0]

    # DT[is.na(new.TRT.gstar), new.TRT.gstar := 0]
    setnames(DT, old = "new.TRT.gstar", new = new.TRT.names[dtheta.i])
  }
  if (!return.allcolumns) {
    DT <- DT[, c(ID, t, new.TRT.names), with = FALSE]
  }
  return(DT)
}

# ----------------------------------------------------------------
#### SUMMARY STATISTICS
# ----------------------------------------------------------------
makeFreqTable <- function(rawFreq){
  ntot <- sum(rawFreq)
  rawFreqPercent <- rawFreq / ntot * 100
  fineFreq <- cbind("Frequency"=rawFreq,"\\%"=round(rawFreqPercent,2),"Cumulative Frequency"=cumsum(rawFreq),"Cumulative \\%"=round(cumsum(rawFreqPercent),2))
  fineFreq <- as.data.frame(fineFreq)
  eval(parse(text=paste('fineFreq <- cbind(',deparse(substitute(rawFreq)),'=factor(c("',paste(names(rawFreq),collapse='","'),'")),fineFreq)',sep="")))
  rownames(fineFreq) <- 1:nrow(fineFreq)
  fineFreq <- as.matrix(fineFreq)
  return(fineFreq)
}

# nocov start
make.table.m0 <- function(S.IPAW, RDscale = "-" , nobs = 0, esti = "IPW", t.period, t.value, se.RDscale.Sdt.K, TMLE = FALSE){
  if (missing(se.RDscale.Sdt.K)) {
    se.RDscale.Sdt.K <- matrix(NA, nrow = length(S.IPAW), ncol = length(S.IPAW))
    colnames(se.RDscale.Sdt.K) <- names(S.IPAW)
    rownames(se.RDscale.Sdt.K) <- names(S.IPAW)
  }
  dtheta <- names(S.IPAW)
  RDtable <- matrix(NA, nrow = (length(dtheta)-1)*2, ncol = length(dtheta)-1)
  # RDtable <- matrix(NA,nrow=2*3,ncol=3)

  ContrastScale <- ifelse(RDscale,"-","/")

  H0val <- ifelse(RDscale,0,1)
  ## Compare mean bootstrap to point estimates - should be similar
  PYK1.IPAW <- allRDtable <- vector("list",length(S.IPAW))
  names(PYK1.IPAW) <- names(allRDtable) <- names(S.IPAW)
  PYK1.IPAW <- lapply(S.IPAW,function(x,y)return(1-x[y]),y=t.period)

  rownames(RDtable) <- rep(rev(dtheta)[-length(dtheta)],each=2)
  colnames(RDtable) <- dtheta[-length(dtheta)]

  for(rule1 in rev(dtheta)[-length(dtheta)]) {
    for(rule2 in dtheta[ 1:(which(dtheta==rule1)-1) ]){
      RD <- mapply(ContrastScale,PYK1.IPAW[[rule1]],PYK1.IPAW[[rule2]])
      (pt <- round(RD,4))
      (CI <- round(RD+c(qnorm(0.025),-qnorm(0.025))*se.RDscale.Sdt.K[rule1,rule2],4))
      (ptCI <- paste(pt," [",CI[1],";",CI[2],"]",sep=""))
      (pval <- paste("SE=",round(se.RDscale.Sdt.K[rule1,rule2],4),", p=", round(2*pnorm( abs((RD-H0val)/se.RDscale.Sdt.K[rule1, rule2]), lower.tail=FALSE ),2) ,sep=""))
      RDtable[ rownames(RDtable)%in%rule1,rule2] <- c(ptCI,pval)
    }
  }
  RDtable[is.na(RDtable)] <- ""
  RDtable <- cbind("d1 (row) | d2 (col)"=rep(rev(dtheta)[-length(dtheta)],each=2) , RDtable)
  fileText <- paste(esti,ifelse(RDscale,"RD","RR"),sep="")
  captionText2 <- ifelse(RDscale,"differences","ratios")
  captionText3 <- ifelse(RDscale,"$-$","$/$")

  if (esti %in% c("IPAW", "IPW")) {
    estimates <- "Stabilized inverse weighting"
  } else if(esti %in% c("IPAWtrunc","IPWtrunc")){
    estimates <- "Stabilized, truncated inverse weighting"
  } else if (esti %in% "crude") {
    estimates <- "Crude"
  } else {
    estimates <- esti
  }

  model <- "MSM"
  if (esti %in% "crude") model <- "model"

  MSM_text <- paste0("The risk contrasts are derived from a logistic ", model, " for the discrete-time hazards fitted based on ", nobs, " observations. ")
  caption <- paste0(estimates, " estimates of the (cumulative) risk ",
                    captionText2,", $d_1$", captionText3,"$d_2$, evaluated at the time indicator value = ", t.value, ". ",
                    ifelse(!TMLE, MSM_text, ""),
                    "Variance estimates are derived based on the influence curve of the estimator.")
  rownames(RDtable) <- NULL
  return(list(RDtable = RDtable, caption = caption))
}
# nocov end




# nocov start
## RD:
getSE_table_d_by_d <- function(S2.IPAW, IC.Var.S.d, nID, t.period.val.idx, getSEs) {
  se.RDscale.Sdt.K <- matrix(NA, nrow = length(S2.IPAW), ncol = length(S2.IPAW))
  colnames(se.RDscale.Sdt.K) <- names(S2.IPAW)
  rownames(se.RDscale.Sdt.K) <- names(S2.IPAW)
  for (d1.idx in seq_along(names(S2.IPAW))) {
    for (d2.idx in seq_along(names(S2.IPAW))) {
      #### GET SE FOR RD(t)=Sd1(t) - Sd2(t)
      if (getSEs) {
        se.RDscale.Sdt.K[d1.idx, d2.idx] <- getSE.RD.d1.minus.d2(nID = nID,
                                                                 IC.S.d1 = IC.Var.S.d[[d1.idx]][["IC.S"]],
                                                                 IC.S.d2 = IC.Var.S.d[[d2.idx]][["IC.S"]])[t.period.val.idx]


      }
    }
  }
  return(se.RDscale.Sdt.K)
}
# nocov end

# nocov start
#' Risk Difference Estimates and SEs for IPW-MSM
#'
#' Produces table(s) with pair-wise risk differences for all regimens that were used for fitting IPW-MSM.
#' The corresponding SEs are evaluated based on the estimated influence curves (IC).
#' @param MSM Object returned by \code{\link{survMSM}}.
#' @param t.periods.RDs Vector of time-points for evaluation of pairwise risk differences.
#' @param getSEs Evaluate the influence curve based RD estimates of standard errors (SEs) along with point estimates?
#' @return A list with RD tables. One table for each time-point in \code{t.periods.RDs}.
#' @seealso \code{\link{survMSM}} for estimation with MSM.
#' @export
get_MSM_RDs <- function(MSM, t.periods.RDs, getSEs = TRUE) {
  RDs.IPAW.tperiods <- vector(mode = "list", length = length(t.periods.RDs))
  periods_idx <- seq_along(MSM$periods)
  names(RDs.IPAW.tperiods) <- "RDs_for_t" %+% t.periods.RDs
  for (t.idx in seq_along(t.periods.RDs)) {
    t.period.val.idx <- periods_idx[MSM$periods %in% t.periods.RDs[t.idx]]
    se.RDscale.Sdt.K <- getSE_table_d_by_d(MSM$St, MSM$IC.Var.S.d, MSM$nID, t.period.val.idx, getSEs)
    RDs.IPAW.tperiods[[t.idx]] <- make.table.m0(MSM$St,
                                                RDscale = TRUE,
                                                t.period = t.period.val.idx,
                                                t.value = t.periods.RDs[t.idx],
                                                nobs = MSM$nobs,
                                                # nobs = nrow(MSM$wts_data),
                                                esti = MSM$est_name,
                                                se.RDscale.Sdt.K = se.RDscale.Sdt.K)
  }
  return(RDs.IPAW.tperiods)
}

#' Risk Difference Estimates and SEs for a list of TMLE outputs
#'
#' Produces table(s) with pair-wise risk differences for all regimens that were used for fitting TMLE.
#' The corresponding SEs are evaluated based on the estimated influence curves (IC).
#' @param TMLE_list A list of objects returned by \code{\link{fit_iterTMLE}}, \code{\link{fit_TMLE}} or \code{\link{fit_GCOMP}}.
#' @param t.periods.RDs Vector of time-points for evaluation of pairwise risk differences.
#' @return A list with RD tables. One table for each time-point in \code{t.periods.RDs}.
#' @seealso \code{\link{survMSM}} for estimation with MSM.
#' @export
get_TMLE_RDs <- function(TMLE_list, t.periods.RDs) {
  rule_names <- lapply(TMLE_list, "[[", "rule_name")
  names(TMLE_list) <- rule_names
  new_TMLE_list <- list()
  new_TMLE_list$St <- lapply(lapply(TMLE_list, "[[", "estimates"), '[[', "St.TMLE")
  new_TMLE_list$IC.Var.S.d <- lapply(TMLE_list, "[[", "IC.Var.S.d")
  new_TMLE_list$periods <- TMLE_list[[1]]$periods
  new_TMLE_list$est_name <- TMLE_list[[1]]$est_name
  new_TMLE_list$wts_data <- rbindlist(lapply(TMLE_list, "[[", "wts_data"))
  new_TMLE_list$nID <- TMLE_list[[1]]$nID

  RDs.TMLE.tperiods <- vector(mode = "list", length = length(t.periods.RDs))
  periods_idx <- seq_along(new_TMLE_list$periods)
  names(RDs.TMLE.tperiods) <- "RDs_for_t" %+% t.periods.RDs
  for (t.idx in seq_along(t.periods.RDs)) {
    t.period.val.idx <- periods_idx[new_TMLE_list$periods %in% t.periods.RDs[t.idx]]
    se.RDscale.Sdt.K <- getSE_table_d_by_d(new_TMLE_list$St, new_TMLE_list$IC.Var.S.d, new_TMLE_list$nID, t.period.val.idx, getSEs = TRUE)
    RDs.TMLE.tperiods[[t.idx]] <- make.table.m0(new_TMLE_list$St,
                                                RDscale = TRUE,
                                                t.period = t.period.val.idx,
                                                t.value = t.periods.RDs[t.idx],
                                                nobs = nrow(new_TMLE_list$wts_data),
                                                esti = new_TMLE_list$est_name,
                                                se.RDscale.Sdt.K = se.RDscale.Sdt.K,
                                                TMLE = TRUE)
  }
  return(RDs.TMLE.tperiods)
}
# nocov end
