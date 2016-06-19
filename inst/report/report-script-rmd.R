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
f_plot_survest <- function(surv_res_est, t_int_sel, y_lab, x_lab, miny, x_legend, y_legend) {
  ptsize <- 0.7
  counter <- 0
  if (missing(y_lab)) y_lab <- ""
  if (missing(x_lab)) x_lab <- "Follow-up period since study entry"
  if (missing(t_int_sel)) t_int_sel <- seq_along(surv_res_est[[1]])
  if (missing(miny)) miny <- min(unlist(lapply(surv_res_est, function(x) min(x[t_int_sel], na.rm = TRUE))))
  if (missing(x_legend)) x_legend <- (max(t_int_sel, na.rm = TRUE) - min(t_int_sel, na.rm = TRUE)) * 2/3 + min(t_int_sel, na.rm = TRUE)
  if (missing(y_legend)) y_legend <- (1 - miny) * 4/5 + miny
  for(d.j in names(surv_res_est)){
    counter <- counter + 1
    plot(t_int_sel, surv_res_est[[d.j]][t_int_sel], col = counter, type = 'b', cex = ptsize, ylim = c(miny, 1), ylab = y_lab, xlab = x_lab)
    par(new=TRUE)
  }
  legend(x_legend, y_legend, legend = names(surv_res_est), col = c(1:length(names(surv_res_est))), cex = ptsize, pch = 1)
}

#' # Model fits for propensity scores
#'
#' Number of unique IDs in the input data:
{{prettyNum(OData$nuniqueIDs, big.mark = ",", scientific = FALSE)}}
#'
#' Number of person-time observations in the input data:
{{prettyNum(OData$nobs, big.mark = ",", scientific = FALSE)}}
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
#' # Distribution of the weights

#+ echo=FALSE
if (!missing(WTtables)) {
  # IPAWdists <- get_wtsummary(MSM$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150))
  # IPAWdistByRule <- IPAWdists$summary.table.byrule
  pander::set.caption("Distribution of the stabilized IPA weights for all rule-person-time observations")
  pander::pander(WTtables$summary.table, justify = c('right', rep("left",ncol(WTtables$summary.table)-1)))
}

#+ echo=FALSE
if (!missing(WTtables) & !is.null(WTtables$summary.DT.byrule)) {
  pander::set.caption("Counts of the stabilized IPA weights by each rule")
  pander::pander(WTtables$summary.DT.byrule, justify = c('right', rep("left",ncol(WTtables$summary.DT.byrule)-1)))
}

#'\pagebreak
#'
#' # MSM fits

#+ echo=FALSE
MSM.fit <- MSM$MSM.fit
output.MSM <- round(MSM.fit$coef,2)
output.MSM <- cbind("Terms" = names(MSM.fit$coef), output.MSM)
colnames(output.MSM) <- c("Terms",ifelse(MSM$trunc_weights == Inf && MSM$use_weights, "IPAW", ifelse(MSM$trunc_weights < Inf && MSM$use_weights, "truncated IPAW", "no weights")))
rownames(output.MSM) <- NULL
pander::set.caption("Coefficients of MSM")
pander::pander(output.MSM, justify = c('right', 'left'))

#'\pagebreak
#'
#' # Survival estimates

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "Survival Estimates.\\label{fig:survPlot}"
sysArg <- list()
sysArg$surv_res_est <- SurvByRegimen
userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
do.call(f_plot_survest, sysArg)

#'\pagebreak
#'
#' # RD tables

#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
if (!missing(RDtables)) {
  for (RDtable in RDtables) {
    pander::set.caption(RDtable$caption)
    pander::pander(RDtable$RDtable)
  }
}

#+ include=FALSE
panderOptions('knitr.auto.asis', TRUE)
