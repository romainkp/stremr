test.iterTMLE.10Kdata <- function() {
  # options(stremr.verbose = TRUE)
  options(stremr.verbose = FALSE)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  # set_all_stremr_options(estimator = "speedglm__glm")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # Odat_DT <- Odat_DT[ID %in% (1:100), ]
  setkeyv(Odat_DT, cols = c("ID", "t"))

  # ---------------------------------------------------------------------------
  # Define some summaries (lags C[t-1], A[t-1], N[t-1])
  # ---------------------------------------------------------------------------
  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- lagnodes %+% ".tminus1"
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

  # ----------------------------------------------------------------
  # IMPORT DATA
  # ----------------------------------------------------------------
  # options(stremr.verbose = TRUE)

  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  # ------------------------------------------------------------------
  # Fit propensity scores for Treatment, Censoring & Monitoring
  # ------------------------------------------------------------------
  gform_TRT <- c("TI ~ CVD + highA1c")
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))

  gform_CENS <- c("C ~ highA1c + t")
  gform_MONITOR <- "N ~ 1"
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  # ---------------------------------------------------------------------------------------------------------
  # Iterative TMLE
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(2)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  # params = list(fit.package = "speedglm", fit.algorithm = "glm")
  # params = list(fit.package = "h2o", fit.algorithm = "RF", ntrees = 100,
  #               learn_rate = 0.05, sample_rate = 0.8,
  #               col_sample_rate = 0.8, balance_classes = TRUE)
  # Qstratify <- c("TI == 0 & CVD == 0", "TI == 1 & CVD == 0", "TI == 0 & CVD == 1", "TI == 1 & CVD == 1")
  # Qstratify = Qstratify,
  # models = params,
  iterTMLE_est1a <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  iterTMLE_est1a$estimates[]
#             est_name t      risk      surv ALLsuccessTMLE nFailedUpdates   type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.1144987 0.8855013          FALSE              5 pooled    0.1123336    0.8876664 1.897894e-05 0.004356482 gTI.dhigh
  # models = params,
  iterTMLE_est1b <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  iterTMLE_est1b$estimates[]
#               est_name t       risk      surv ALLsuccessTMLE nFailedUpdates   type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.04205589 0.9579441          FALSE              5 pooled    0.0541523    0.9458477 3.802643e-05 0.006166557  gTI.dlow
  # models = params,
  iterTMLE_est2a <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  iterTMLE_est2a$estimates[]
#             est_name t      risk      surv ALLsuccessTMLE nFailedUpdates       type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.1123625 0.8876375          FALSE              5 stratified     0.112169     0.887831 1.897208e-05 0.004355695 gTI.dhigh
  # models = params,
  iterTMLE_est2a_nodataadapt <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE, adapt_stop = FALSE)
  iterTMLE_est2a_nodataadapt$estimates[]
#             est_name t      risk      surv ALLsuccessTMLE nFailedUpdates       type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.1123625 0.8876375          FALSE              5 stratified    0.1118552    0.8881448 1.897197e-05 0.004355682 gTI.dhigh
  # models = params,
  iterTMLE_est2b <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  iterTMLE_est2b$estimates[]
#             est_name t       risk      surv ALLsuccessTMLE nFailedUpdates       type iterTMLErisk iterTMLEsurv    TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.04856046 0.9514395          FALSE              5 stratified   0.05410295    0.9458971 3.80063e-05 0.006164925  gTI.dlow

  # -------------------------------------------------------
  # Changing the factor in the adaptive stopping criteria (multiplier of 1/ (f*sqrt(N)))
  # -------------------------------------------------------
  # The higher the factor -> the more iterations are needed:
  # models = params,
  iterTMLE_est2a_fact50 <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE, adapt_stop_factor = 50)
  iterTMLE_est2a_fact50$estimates[]
#             est_name t      risk      surv ALLsuccessTMLE nFailedUpdates       type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.1123625 0.8876375          FALSE              5 stratified    0.1118552    0.8881448 1.897197e-05 0.004355682 gTI.dhigh

  # the lower the factor -> the fewer iterations are needed:
  # models = params,
  iterTMLE_est2b_fact2 <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE, adapt_stop_factor = 2)
  iterTMLE_est2b_fact2$estimates[]
#             est_name t       risk      surv ALLsuccessTMLE nFailedUpdates       type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE
# 1: GCOMP & Iter.TMLE 4 0.04856046 0.9514395          FALSE              5 stratified   0.05043029    0.9495697 3.802771e-05 0.006166662

  # -------------------------------------------------------
  # NON-DATA-ADAPTIVE stopping criteria
  # -------------------------------------------------------
  # models = params,
  iterTMLE_est2c <- fit_iterTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE, adapt_stop = FALSE)
  iterTMLE_est2c$estimates[]
#             est_name t       risk      surv ALLsuccessTMLE nFailedUpdates       type iterTMLErisk iterTMLEsurv     TMLE_Var     TMLE_SE rule.name
# 1: GCOMP & Iter.TMLE 4 0.04856046 0.9514395          FALSE              5 stratified   0.05478108    0.9452189 3.800412e-05 0.006164748  gTI.dlow
}

