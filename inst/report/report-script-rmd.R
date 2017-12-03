#' ---
#' title: "`r title`"
#' author: "`r author`"
#' date: "`r Sys.Date()`"
#' ---

#+ setup, include=FALSE
require("knitr")
require("pander")
# require("gridisl")
opts_chunk$set(fig.path = figure.dir)
panderOptions("table.split.table", Inf)

print_model_info <- function(model_summary, model_stack) {
  cat("\n\n"); cat("###"); cat("Model Performance"); cat("\n\n");

  pander::pander(model_summary)
  MSEtab <- model_stack$getMSEtab
  try(pander::pander(MSEtab, caption = "Overall Performance by Model"))
  
  if (is(model_stack, "Lrnr_base")) {
    fit <- model_stack$fit_object
    pander::pander(fit)
  } else {
    grids <- model_stack$get_modelfits_grid()
    for (grid in grids) {
      if (is.data.frame(grid) || is.data.table(grid)) {
        grid <- grid[ , names(grid)[!(names(grid) %in% c("glob_params", "xgb_fit", "fit", "params"))], with = FALSE]
        try(pander::pander(grid, caption = "XGB Grid"))
        # cat(paste(capture.output(print(grid)), collapse = '\n\n'))
      } else {
        # cat(paste(capture.output(print(grid))[-2], collapse = '\n\n'))
        try(pander::pander(grid))
        # try(pander::pander(paste(capture.output(print(grid))[-2], collapse = '\n')))
      }
    }    
  }
  # cat("\n\n"); cat("###"); cat("Best Model"); cat("\n\n");
  # # best_model <- model_stack$get_best_models(K=1)[[1]]
  # best_model <- model_stack$get_overall_best_model()[[1]]
  # gridisl::print_tables(best_model)
}

#'
#' Number of unique independent units in the input data:
{{prettyNum(nuniqueIDs, big.mark = ",", scientific = FALSE)}}
#'
#' Number of person-time observations in the input data:
{{prettyNum(nobs, big.mark = ",", scientific = FALSE)}}
#'
#' Total number of unique time-points in the input data:
{{prettyNum(nuniquets, big.mark = ",", scientific = FALSE)}}
#'
#' # Model fits for propensity scores
#'
#' ## Model(s) for censoring variable(s):

#+ echo=FALSE, warning=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
set.alignment('left', row.names = 'right')
if (!skip.modelfits) {
  for (reg.model.idx in seq_along(model_fits_gC)) {
    model_summary <- model_summaries_gC[[reg.model.idx]]
    model_stack <- model_fits_gC[[reg.model.idx]]
    print_model_info(model_summary, model_stack)
  }
}

#' ## Model(s) for exposure variable(s):

#+ echo=FALSE, warning=FALSE, results='asis'
if (!skip.modelfits) {
  for (reg.model.idx in seq_along(model_fits_gA)) {
    # cat("\n\n"); cat("###"); cat("Model Summary"); cat("\n\n");
    # pander::pander(model_summaries_gA[[reg.model.idx]])
    # best_model <- model_fits_gA[[reg.model.idx]]$get_best_models(K=1)[[1]]
    # gridisl::print_tables(best_model)
    # # print(reg.model$get_best_models(), only.coefs = only.coefs)
    model_summary <- model_summaries_gA[[reg.model.idx]]
    model_stack <- model_fits_gA[[reg.model.idx]]
    print_model_info(model_summary, model_stack)
  }
}

#' ## Model(s) for monitoring variable(s):

#+ echo=FALSE, warning=FALSE, results='asis'
if (!skip.modelfits) {
  for (reg.model.idx in seq_along(model_fits_gN)) {
    # cat("\n\n"); cat("###"); cat("Model Summary"); cat("\n\n");
    # pander::pander(model_summaries_gN[[reg.model.idx]])
    # best_model <- model_fits_gN[[reg.model.idx]]$get_best_models(K=1)[[1]]
    # gridisl::print_tables(best_model)
    # # print(reg.model$get_best_models(), only.coefs = only.coefs)
    model_summary <- model_summaries_gN[[reg.model.idx]]
    model_stack <- model_fits_gN[[reg.model.idx]]
    print_model_info(model_summary, model_stack)
  }
}

#+ include=FALSE
panderOptions('knitr.auto.asis', TRUE)

#'\pagebreak
#'
#' `r ifelse(!missing(WTtables),'# Distribution of the weights','')`

#+ echo=FALSE
if (!missing(WTtables)) {
  if (!is.null(WTtables)) {
    pander::set.caption("Distribution of the stabilized IPA weights for all rule-person-time observations")
    pander::pander(WTtables$summary.table, justify = c('right', rep("left",ncol(WTtables$summary.table)-1)))
  }
}

#+ echo=FALSE
if (!missing(WTtables)) {
  if (!is.null(WTtables) && !is.null(WTtables$summary.DT.byrule)) {
    pander::set.caption("Counts of the stabilized IPA weights by each rule")
    pander::pander(WTtables$summary.DT.byrule, justify = c('right', rep("left",ncol(WTtables$summary.DT.byrule)-1)))
  }
}

#'\pagebreak
#'
#' `r ifelse(!missing(FUPtables),'# Distribution of the follow-up times','')`

#+ echo=FALSE
if (!missing(FUPtables)) {
  if (!is.null(FUPtables)) {
    rules <- unique(FUPtables[["rule.name"]])
    for (T.rule in rules) {
      one_ruleID <- FUPtables[(rule.name %in% eval(T.rule)), max.t]
      hist(one_ruleID, main = "Maximum follow-up period for TRT/MONITOR rule: " %+% T.rule)
    }
  }
}

#+ echo=FALSE, results='asis'
if (!missing(FUPtables)) {
  if (!is.null(FUPtables)) {
    rules <- unique(FUPtables[["rule.name"]])
    for (T.rule in rules) {
      one_ruleID <- FUPtables[(rule.name %in% eval(T.rule)), max.t]
      panderOptions('knitr.auto.asis', FALSE)
      followupTimes <- table(one_ruleID)
      followupTimes <- makeFreqTable(followupTimes)
      # followupTimes <- makeFreqTable(table(one_ruleID))
      pander::pander(followupTimes, caption = "Distribution of the total follow-up time for TRT/MONITOR rule: " %+% T.rule)
      pander::pander(summary(one_ruleID), caption = "Min/Max/Quantiles for the total follow-up time for TRT/MONITOR rule: " %+% T.rule)
      panderOptions('knitr.auto.asis', TRUE)
    }
  }
}

#'\pagebreak
#'
#' `r ifelse(!missing(NPMSM) && plotKM,'# Survival with Kaplan-Meier','')`

#+ echo=FALSE, warning = FALSE, fig.width=8, fig.height=5, fig.cap = "Survival with KM.\\label{fig:survPlotGCOMP}"
if (!missing(NPMSM) && plotKM) {
  est_name <- "KM"
  est_obj <- NPMSM

  sysArg <- list()
  if ("estimates" %in% names(est_obj)) est_obj <- list(res = est_obj)
  est_obj <- lapply(est_obj, '[[', "estimates")

  if (!use_ggplot) {

    sysArg <- f_obtain_St(sysArg, est_obj, optArgReport, est_name = "St."%+%est_name, t_name = "time")
    do.call(f_plot_survest, sysArg)

  } else {

    sysArg[["estimates"]] <- est_obj
    sysArg[["surv_name"]] <- "St."%+%est_name

    userArg <- intersect(names(formals(ggsurv)), names(optArgReport)) # captures optional arguments given by user for customizing report
    if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
    do.call(ggsurv, sysArg)

  }
}

#'\pagebreak
#'
#' `r ifelse(!missing(NPMSM),'# Survival with IPW-Adjusted Kaplan-Meier (Non-Parametric / Saturated MSM for Hazard)','')`

#+ echo=FALSE, warning=FALSE, fig.width=8, fig.height=5, fig.cap = "Survival with IPW-Adjusted KM.\\label{fig:survPlotGCOMP}"
if (!missing(NPMSM)) {
  est_name <- "NPMSM"
  est_obj <- NPMSM

  sysArg <- list()
  if ("estimates" %in% names(est_obj)) est_obj <- list(res = est_obj)
  est_obj <- lapply(est_obj, '[[', "estimates")

  if (!use_ggplot) {

    sysArg <- f_obtain_St(sysArg,est_obj, optArgReport, est_name = "St."%+%est_name, t_name = "time")
    do.call(f_plot_survest, sysArg)

  } else {

    sysArg[["estimates"]] <- est_obj
    userArg <- intersect(names(formals(ggsurv)), names(optArgReport)) # captures optional arguments given by user for customizing report
    if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
    do.call(ggsurv, sysArg)

  }

}

#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
if (!missing(NPMSM) && printEstimateTables) {
  for (NPMSMtab in surv_tables) {
    pander::set.caption("NPMSM results for rule '" %+% NPMSMtab[["rule.name"]][1] %+% "'")
    pander::pander(data.frame(NPMSMtab))
  }
}
panderOptions('knitr.auto.asis', TRUE)

#'\pagebreak
#'
#' `r ifelse(!missing(MSM),'# IPW-based MSM fits','')`

#+ echo=FALSE
if (!missing(MSM)) {
  MSM.fit <- MSM$MSM.fit
  output.MSM <- round(MSM.fit$coef,2)
  output.MSM <- cbind("Terms" = names(MSM.fit$coef), output.MSM)
  # colnames(output.MSM) <- c("Terms",ifelse(MSM$trunc_weights == Inf && MSM$use_weights, "IPAW", ifelse(MSM$trunc_weights < Inf && MSM$use_weights, "truncated IPAW", "no weights")))
  colnames(output.MSM) <- c("Terms", "Estimates")
  rownames(output.MSM) <- NULL
  pander::set.caption("Coefficients of logistic regression model for hazards")
  pander::pander(output.MSM, justify = c('right', 'left'))
}

#'\pagebreak
#'
#' `r ifelse(!missing(MSM),'# Survival with IPW-based Marginal Structural Model (MSM)','')`

#+ echo=FALSE, warning = FALSE, fig.width=8, fig.height=5, fig.cap = "IPW-MSM Survival.\\label{fig:survPlotIPW}"
if (!missing(MSM)) {
  est_name <- "MSM"
  est_obj <- MSM

  sysArg <- list()

  if (!use_ggplot) {

    sysArg$surv_list <- lapply(est_obj[["estimates"]], "[[", "St."%+%est_name)
    sysArg$t <- est_obj[["periods"]]
    userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
    if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
    do.call(f_plot_survest, sysArg)

  } else {

    sysArg[["estimates"]] <- est_obj[["estimates"]]
    userArg <- intersect(names(formals(ggsurv)), names(optArgReport)) # captures optional arguments given by user for customizing report
    if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
    do.call(ggsurv, sysArg)

  }

}

#'\pagebreak
#'
#' `r ifelse(!missing(MSM.RDtables),'# IPW-MSM RD Tables','')`

#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
if (!missing(MSM.RDtables)) {
  for (RDtable in MSM.RDtables) {
    pander::set.caption(RDtable$caption)
    pander::pander(RDtable$RDtable)
  }
}
panderOptions('knitr.auto.asis', TRUE)

#' `r ifelse(!missing(GCOMP),'# Survival with Sequential G-Computation','')`

#+ echo=FALSE, warning = FALSE, fig.width=8, fig.height=5, fig.cap = "G-Computation Survival.\\label{fig:survPlotGCOMP}"
if (!missing(GCOMP)) {
  est_name <- "GCOMP"
  est_obj <- GCOMP

  sysArg <- list()
  if ("estimates" %in% names(est_obj)) est_obj <- list(res = est_obj)
  est_obj <- lapply(est_obj, '[[', "estimates")

  if (!use_ggplot) {

    sysArg <- f_obtain_St(sysArg,est_obj, optArgReport, est_name = "St."%+%est_name, t_name = "time")
    do.call(f_plot_survest, sysArg)

  } else {

    sysArg[["estimates"]] <- est_obj
    userArg <- intersect(names(formals(ggsurv)), names(optArgReport)) # captures optional arguments given by user for customizing report
    if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
    do.call(ggsurv, sysArg)

  }
}

#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
if (!missing(GCOMP) && printEstimateTables) {
  GCOMP.St <- lapply(GCOMP, '[[', "estimates")
  for (GCOMPtab in GCOMP.St) {
    pander::set.caption("GCOMP results for rule '" %+% GCOMPtab[["rule.name"]][1] %+% "'")
    pander::pander(data.frame(GCOMPtab))
  }
}
panderOptions('knitr.auto.asis', TRUE)

#' `r ifelse(!missing(TMLE),'# Survival with Targeted Maximum Likelihood (TMLE)','')`

#+ echo=FALSE, warning = FALSE, fig.width=8, fig.height=5, fig.cap = "TMLE Survival.\\label{fig:survPlotTMLE}"
if (!missing(TMLE)) {
  # sysArg <- f_obtain_St(sysArg,TMLE, optArgReport, est_name = "St.TMLE", t_name = "time")
  # do.call(f_plot_survest, sysArg)
  est_name <- "TMLE"
  est_obj <- TMLE

  sysArg <- list()
  if ("estimates" %in% names(est_obj)) est_obj <- list(res = est_obj)
  est_obj <- lapply(est_obj, '[[', "estimates")

  if (!use_ggplot) {

    sysArg <- f_obtain_St(sysArg,est_obj, optArgReport, est_name = "St."%+%est_name, t_name = "time")
    do.call(f_plot_survest, sysArg)

  } else {

    sysArg[["estimates"]] <- est_obj
    userArg <- intersect(names(formals(ggsurv)), names(optArgReport)) # captures optional arguments given by user for customizing report
    if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
    do.call(ggsurv, sysArg)

  }
}

#'\pagebreak
#'
#' `r ifelse(!missing(TMLE.RDtables),'# TMLE RD Tables','')`

#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
if (!missing(TMLE.RDtables)) {
  for (RDtable in TMLE.RDtables) {
    pander::set.caption(RDtable$caption)
    pander::pander(RDtable$RDtable)
  }
}
panderOptions('knitr.auto.asis', TRUE)

#'\pagebreak
#'
#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
if (!missing(TMLE) && printEstimateTables) {
  TMLE.St <- lapply(TMLE, '[[', "estimates")
  for (TMLEtab in TMLE.St) {
    pander::set.caption("TMLE results for rule '" %+% TMLEtab[["rule.name"]][1] %+% "'")
    pander::pander(data.frame(TMLEtab))
  }
}
panderOptions('knitr.auto.asis', TRUE)
