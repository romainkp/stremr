# ------------------------------------------------------------------------------------------------------
# Testing stremr building blocks with standard glm.fit
# ------------------------------------------------------------------------------------------------------
test.buildingblocks <- function() {
  require("data.table")
  require("stremr")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  # set_all_stremr_options(estimator = "glm__glm")

  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  # Nsize <- 1000
  # OdataNoCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  data(OdataNoCENS)
  OdataDT <- as.data.table(OdataNoCENS, key=c("ID", "t"))

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
               estimates

  res.test.AMKE.IPTW.1 <- c(0.9564462, 0.9403990, 0.9250282, 0.9250282, 0.9010454, 0.8873632, 0.8850073,
                            0.8678221, 0.8488353, 0.8401469, 0.8401469, 0.8304123, 0.7930334, 0.7752481,
                            0.7752481, 0.7340285, 0.7340285)
  checkTrue(all(abs(res.test.AMKE.IPTW.1 - AKME.St.1[["St.NPMSM"]]) < (10^-5)))

  res.test.AMKE.KM.1 <- c(0.9526627, 0.9349112, 0.9230769, 0.9230769, 0.8994083, 0.8757396, 0.8639053,
                          0.8520710, 0.8224852, 0.8165680, 0.8165680, 0.8047337, 0.7692308, 0.7573964,
                          0.7573964, 0.7218935, 0.7218935)
  checkTrue(all(abs(res.test.AMKE.KM.1 - AKME.St.1[["St.KM"]]) < (10^-5)))
  # ----------------------------------------------------------------------
  # Bounded IPW
  # ----------------------------------------------------------------------
  IPW.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
               directIPW(OData)
  # IPW.St.1[]
  res.test.IPW.St.1 <- c(0.9564462, 0.9519073, 0.9497037, 0.9643112, 0.9498030, 0.9502507, 0.9631751,
                         0.9559085, 0.9487348, 0.9553832, 0.9698508, 0.9683569, 0.9347756, 0.9341488,
                         0.9552954, 0.9183222)
  checkTrue(all(abs(res.test.IPW.St.1 - IPW.St.1[["S.t.n"]]) < (10^-5)))

  # ----------------------------------------------------------------------
  # IPW-MSM for hazard
  # ----------------------------------------------------------------------
  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
  wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
  survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, tbreaks = c(1:8,12,16)-1,)
  # survMSM_res$St
  res.test <- list()
  res.test$TI0 <- c(0.9957857, 0.9827720, 0.9653188, 0.9653188, 0.9653188, 0.9653188, 0.9641698,
                    0.9625935, 0.9625935, 0.9625935, 0.9625935, 0.9625935, 0.9625935, 0.9625935,
                    0.9625935, 0.9625935)
  res.test$TI1 <- c(0.9564462, 0.9403990, 0.9250282, 0.9250277, 0.9010449, 0.8873627, 0.8850068,
                    0.8678217, 0.8595666, 0.8513901, 0.8432913, 0.8352696, 0.8085728, 0.7827293,
                    0.7577119, 0.7334940)
  # res <- all.equal(res.test, survMSM_res$St)
  checkTrue(all(abs(survMSM_res$St$TI0 - res.test$TI0) < (10^-5)))
  checkTrue(all(abs(survMSM_res$St$TI1 - res.test$TI1) < (10^-5)))

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
  survNPIPW_ests <- survNP_res$estimates
  survNPIPW_ests[]

  survDirIPW_ests <- directIPW(wts.DT.1, OData)
  survDirIPW_ests[]

  survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, tbreaks = c(1:8,12,16)-1,)
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
  t.surv <- c(0:2)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est[]

  # stratified modeling by rule followers only:
  tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est[]
}
