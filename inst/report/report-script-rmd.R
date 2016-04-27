#' ---
#' title: "stremr Analysis Report"
#' author: "Insert Author"
#' date: "`r Sys.Date()`"
#' ---

#+ setup, include=FALSE
library("knitr")
opts_chunk$set(fig.path = 'figure/stremr-')


#+ echo=FALSE
f_plot_survest <- function(surv_res_est, t_int_sel, y_lab, miny) {
  ptsize <- 0.7
  counter <- 0
  if (missing(y_lab)) y_lab <- ""
  if (missing(t_int_sel)) t_int_sel <- seq_along(surv_res_est[[1]])
  if (missing(miny)) miny <- min(unlist(lapply(surv_res_est, min)))

  for(d.j in names(surv_res_est)){
    counter <- counter+1
    plot(t_int_sel,surv_res_est[[d.j]][t_int_sel],col=counter,type='b',cex=ptsize,ylim=c(miny,1),
      ylab = y_lab, xlab="Quarter since study entry")
    par(new=TRUE)
  }
  legend(12,0.96,legend=names(surv_res_est),col=c(1:length(names(surv_res_est))), cex=ptsize, pch=1)
}


#' # Model fits for propensity scores

#' ### Model(s) for censoring variable(s)
#+ echo=FALSE, results='asis'
panderOptions('knitr.auto.asis', FALSE)
for (reg.model in fitted.coefs.gC) {
  print(pander::set.caption("Regression: " %+% reg.model$regression))
  print(pander::pander(reg.model$coef, justify = c('right', 'center')))
}

#' ### Model(s) for exposure variable(s)
#+ echo=FALSE, results='asis'
# pander::set.caption("Regression: " %+% fitted.coefs.gA$regression)
# pander::pander(fitted.coefs.gA$coef, justify = c('right', 'center'))
for (reg.model in fitted.coefs.gA) {
  print(pander::set.caption("Regression: " %+% reg.model$regression))
  print(pander::pander(reg.model$coef, justify = c('right', 'center')))
}


#' ### Model(s) for monitoring variable(s)
#+ echo=FALSE, results='asis'
# pander::set.caption("Regression: " %+% fitted.coefs.gN$regression)
# pander::pander(fitted.coefs.gN$coef, justify = c('right', 'center'))
for (reg.model in fitted.coefs.gN) {
  print(pander::set.caption("Regression: " %+% reg.model$regression))
  print(pander::pander(reg.model$coef, justify = c('right', 'center')))
}

#+ include=FALSE
panderOptions('knitr.auto.asis', TRUE)

#' # Distribution of the Weights

#+ echo=FALSE
MSM$IPAWdist

#+ echo=FALSE
pander::set.caption("Distribution of the weights")
pander::pander(MSM$IPAWdist, justify = c('right', rep("center",ncol(MSM$IPAWdist)-1)))

#' # MSM fits
#+ echo=FALSE
pander::set.caption("Coefficients of MSM")
pander::pander(MSM$output.MSM, justify = c('right', 'center'))


#' # Survival estimates

#+ echo=FALSE, fig.width=5, fig.height=5, fig.cap = "Survival Estimates.\\label{fig:survPlot}"
f_plot_survest(Surv.byregimen)