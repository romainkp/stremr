# ------------------------------------------------------------------------------------------------------
# Testing stremr building blocks with standard glm.fit
# ------------------------------------------------------------------------------------------------------
test.buildingblocks <- function() {
  options(stremr.verbose = TRUE)
  require("data.table")
  set_all_stremr_options(fit.package = "glm", fit.algorithm = "glm")

  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  # Nsize <- 1000
  # OdataNoCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  data(OdataNoCENS)
  OdataDT <- as.data.table(OdataNoCENS, key=c(ID, t))

  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataDT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]

  # Define intervention (always treated):
  OdataDT[, ("TI.set1") := 1L]
  OdataDT[, ("TI.set0") := 0L]

  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform_CENS <- "C ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  stratify_CENS <- list(C=c("t < 16", "t == 16"))
  OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "N.tminus1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  checkException(wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1"))

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR, stratify_CENS = stratify_CENS)

  # ----------------------------------------------------------------------
  # IPW Ajusted KM or Saturated MSM
  # ----------------------------------------------------------------------
  require("magrittr")
  AKME.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
               survNPMSM(OData) %$%
               IPW_estimates
  AKME.St.1

  # ----------------------------------------------------------------------
  # Bounded IPW
  # ----------------------------------------------------------------------
  IPW.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
               survDirectIPW(OData)
  IPW.St.1[]

  # ----------------------------------------------------------------------
  # IPW-MSM for hazard
  # ----------------------------------------------------------------------
  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
  wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
  survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, t_breaks = c(1:8,12,16)-1,)
  survMSM_res$St

  # --------------------------------
  # EXAMPLE 2:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
  wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
  survNP_res <- survNPMSM(wts.DT.1, OData)
  survNPIPW_ests <- survNP_res$IPW_estimates
  survNPIPW_ests[]

  survDirIPW_ests <- survDirectIPW(wts.DT.1, OData)
  survDirIPW_ests[]

  survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, t_breaks = c(1:8,12,16)-1,)
  names(survMSM_res)
  survMSM_res$St

  # --------------------------------
  # Test Optional weights:
  # --------------------------------
  checkException(
    survNP_res_addedWts <- survNPMSM(wts.DT.1, OData, weights = 1)
    )
  addedWts_DT <- OdataDT[, c("ID", "t"), with = FALSE]
  addedWts_DT[, new.wts := 1L]
  survNP_res_addedWts <- survNPMSM(wts.DT.1, OData, weights = addedWts_DT)

  # --------------------------------
  # Sequential G-COMP:
  # --------------------------------
  t.surv <- c(0:4)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est[]

  # stratified modeling by rule followers only:
  tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est[]

  # --------------------------------
  res <- stremr(OdataDT, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"), intervened_TRT = "TI.set1",
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
          # noCENScat = 0L)
  res$IPW_estimates

  checkTrue(all.equal(res$IPW_estimates[["St.IPTW"]], survNPIPW_ests[["St.IPTW"]]))
  checkTrue(all.equal(res$IPW_estimates[["St.KM"]], survNPIPW_ests[["St.KM"]]))
  # res$dataDT

  # --------------------------------
  # EXAMPLE 2:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify_CENS <- rep(list(strat.str), 3)
  names(stratify_CENS) <- c("C", "TI", "N")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  wts.DT <- getIPWeights(OData = OData) # , gstar_TRT = ..., gstar_MONITOR = ...)
  survNP_res <- survNPMSM(wts.DT, OData)
  survNPIPW_ests <- survNP_res$IPW_estimates

  res <- stremr(OdataDT, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
          gform_TRT = gform_TRT,
          gform_MONITOR = gform_MONITOR)
          # noCENScat = 0L)
  # res$IPW_estimates
  # res$dataDT
  checkTrue(all.equal(res$IPW_estimates[["St.IPTW"]], survNPIPW_ests[["St.IPTW"]]))
  checkTrue(all.equal(res$IPW_estimates[["St.KM"]], survNPIPW_ests[["St.KM"]]))

  options(stremr.verbose = FALSE)
}
