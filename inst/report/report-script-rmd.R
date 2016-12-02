#' ---
#' title: "`r title`"
#' author: "`r author`"
#' date: "`r Sys.Date()`"
#' ---

#+ setup, include=FALSE
require("knitr")
require("pander")
opts_chunk$set(fig.path = figure.dir)
panderOptions("table.split.table", Inf)

#+ echo=FALSE, include=FALSE
f_plot_survest <- function(surv_list, t, t_int_sel, y_lab, x_lab, miny, x_legend, y_legend, cex = 0.7) {
  # ptsize <- 0.7
  # ptsize <- 0.4
  counter <- 0
  if (missing(y_lab)) y_lab <- ""
  if (missing(x_lab)) x_lab <- "Follow-up period since study entry"
  if (missing(t)) t <- seq_along(surv_list[[1]])
  if (missing(t_int_sel)) t_int_sel <- seq_along(t)
  if (missing(miny)) miny <- min(unlist(lapply(surv_list, function(x) min(x[t_int_sel], na.rm = TRUE))))
  # if (missing(x_legend)) x_legend <- (max(t_int_sel, na.rm = TRUE) - min(t_int_sel, na.rm = TRUE)) * 2/3 + min(t_int_sel, na.rm = TRUE)
  # if (missing(y_legend)) y_legend <- (1 - miny) * 4/5 + miny
  for(d.j in names(surv_list)){
    counter <- counter + 1
    plot(as.integer(t[t_int_sel]), surv_list[[d.j]][t_int_sel], col = counter, type = 'b', cex = cex, ylim = c(miny, 1), ylab = y_lab, xlab = x_lab)
    par(new=TRUE)
  }
  if (missing(y_legend)) {
    if (missing(x_legend)) {
      x_legend <- "bottomleft"
    } else {
      if (!is.character(x_legend)) stop("x_legend must be a character when y_legend is unspecified")
    }
    legend(x_legend, legend = names(surv_list), col = c(1:length(names(surv_list))), cex = cex, pch = 1)
  } else {
    legend(x_legend, y_legend, legend = names(surv_list), col = c(1:length(names(surv_list))), cex = cex, pch = 1)
  }
}
f_obtain_TMLE_St <- function(TMLE, optArgReport) {
  sysArg <- list()
  if (is.data.table(TMLE)) TMLE <- list(TMLE_res = TMLE)
  TMLE <- lapply(TMLE, '[[', "estimates")
  sysArg$surv_list <- lapply(TMLE, '[[', 'surv')
  rule.names <- unlist(lapply(TMLE, function(tmle_res) tmle_res[['rule.name']][1]))
  names(sysArg$surv_list) <- rule.names
  sysArg$t <- TMLE[[1]][["t"]]
  userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
  if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
  return(sysArg)
}

#'
#' Number of unique independent units in the input data:
{{prettyNum(nuniqueIDs, big.mark = ",", scientific = FALSE)}}
#'
#' Number of person-time observations in the input data:
{{prettyNum(nobs, big.mark = ",", scientific = FALSE)}}
#'
#' # Model fits for propensity scores
#'
#' ## Model(s) for censoring variable(s):

#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
set.alignment('left', row.names = 'right')
if (!skip.modelfits) {
  for (reg.model in fitted.coefs.gC) {
    print(reg.model, only.coefs = only.coefs)
  }
}

#' ## Model(s) for exposure variable(s):

#+ echo=FALSE, results='asis'
if (!skip.modelfits) {
  for (reg.model in fitted.coefs.gA) {
    print(reg.model, only.coefs = only.coefs)
  }
}

#' ## Model(s) for monitoring variable(s):

#+ echo=FALSE, results='asis'
if (!skip.modelfits) {
  for (reg.model in fitted.coefs.gN) {
    print(reg.model, only.coefs = only.coefs)
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
#' `r ifelse(!missing(NPMSM),'# Survival with IPW-Adjusted Kaplan-Meier (Non-Parametric / Saturated MSM for Hazard)','')`

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "Survival with IPW-Adjusted KM.\\label{fig:survPlotGCOMP}"
if (!missing(NPMSM)) {
  sysArg <- list()
  if (is.data.table(NPMSM)) NPMSM <- list(NPMSM_res = NPMSM)
  surv_tables <- lapply(NPMSM, '[[', 'IPW_estimates')
  sysArg$surv_list <- lapply(surv_tables, '[[', 'St.IPTW')
  rule.names <- unlist(lapply(surv_tables, function(NPMSM_res) NPMSM_res[['rule.name']][1]))
  names(sysArg$surv_list) <- paste0("St.IPTW: ", rule.names)
  sysArg$t <- NPMSM[[1]][["t"]]
  userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
  if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
  do.call(f_plot_survest, sysArg)
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
#' `r ifelse(!missing(NPMSM) && plotKM,'# Survival with Kaplan-Meier','')`

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "Survival with KM.\\label{fig:survPlotGCOMP}"
if (!missing(NPMSM) && plotKM) {
  sysArg <- list()
  if (is.data.table(NPMSM)) NPMSM <- list(NPMSM_res = NPMSM)
  surv_tables <- lapply(NPMSM, '[[', 'IPW_estimates')
  sysArg$surv_list <- lapply(surv_tables, '[[', 'St.KM')
  rule.names <- unlist(lapply(surv_tables, function(NPMSM_res) NPMSM_res[['rule.name']][1]))
  names(sysArg$surv_list) <- paste0("St.KM: ", rule.names)
  sysArg$t <- NPMSM[[1]][["t"]]
  userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
  if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
  do.call(f_plot_survest, sysArg)
}

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

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "IPW-MSM Survival.\\label{fig:survPlotIPW}"
if (!missing(MSM)) {
  sysArg <- list()
  sysArg$surv_list <- MSM$St
  sysArg$t <- MSM$periods
  userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
  if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
  do.call(f_plot_survest, sysArg)
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

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "G-Computation Survival.\\label{fig:survPlotGCOMP}"
if (!missing(GCOMP)) {
  sysArg <- f_obtain_TMLE_St(GCOMP, optArgReport)
  do.call(f_plot_survest, sysArg)
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

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "TMLE Survival.\\label{fig:survPlotTMLE}"
if (!missing(TMLE)) {
  sysArg <- f_obtain_TMLE_St(TMLE, optArgReport)
  do.call(f_plot_survest, sysArg)
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
