# ---------------------------------------------------------------------------
## Test NDE argument "useonly_t_MONITOR":
## -- Will not intervene on g nodes that satisfy the logical expression provided to this arg
# ---------------------------------------------------------------------------
test.speedglm.stochastic.TMLE.NDE.1Kdata <- function() {
  require("stremr")
  options(stremr.verbose = FALSE)
  options(width = 100)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  Odat_DT <- Odat_DT[ID %in% (1:100), ]
  # define intervention on N as 0101010101...
  Odat_DT[, ("N.star.0101") := t%%2]
  setkeyv(Odat_DT, cols = c("ID", "t"))

  ## ---------------------------------------------------------------------------
  ## Define some summaries (lags C[t-1], A[t-1], N[t-1])
  ## ---------------------------------------------------------------------------
  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- lagnodes %+% ".tminus1"
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
  ## ----------------------------------------------------------------
  ## IMPORT DATA
  ## ----------------------------------------------------------------
  set_all_stremr_options(estimator = "speedglm__glm")
  OData <- importData(Odat_DT,
                      ID = "ID", t = "t",
                      covars = c("highA1c", "lastNat1", "lastNat1.factor"),
                      CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome,
                      remove_extra_rows = FALSE)
  ## ------------------------------------------------------------------
  ## Fit propensity scores for Treatment, Censoring & Monitoring
  ## ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))
  gform_CENS <- c("C ~ highA1c + t")
  gform_MONITOR <- "N ~ 1"

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  ## ---------------------------------------------------------------------------------------------------------
  ## IPW-KM with stochastic intervention on MONITOR
  ## ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly")
  surv1.stoch <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly")
  surv2.stoch <- survNPMSM(wts.St.dhigh, OData)

  ## ---------------------------------------------------------------------------------------------------------
  ## IPW-KM with static intervention on MONITOR
  ## ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101")
  surv1.stat <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101")
  surv2.stat <- survNPMSM(wts.St.dhigh, OData)

  ## ---------------------------------------------------------------------------------------------------------
  ## IPW-KM with static intervention on MONITOR under NDE assumption
  ## -- Will not intervene on g nodes that satisfy the logical expression provided to this arg
  ## ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv1.statNDE <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv2.statNDE <- survNPMSM(wts.St.dhigh, OData)

  ## ---------------------------------------------------------------------------------------------------------
  ## TMLE / GCOMP with a stochastic intervention on MONITOR
  ## ---------------------------------------------------------------------------------------------------------
  t.surv <- c(0:2)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]

  # stratified modeling by rule followers only:
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]

  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]

  ## ---------------------------------------------------------------------------------------------------------
  ## TMLE / GCOMP with a static intervention on MONITOR under NDE assumption
  ## ---------------------------------------------------------------------------------------------------------
  t.surv <- c(0:2)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                            useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]
  # stratified modeling by rule followers only:
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]
}
