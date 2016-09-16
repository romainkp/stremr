# require("R6")

test.buildingblocks <- function() {
  require("data.table")
  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  # Nsize <- 1000
  # OdataNoCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  data(OdataNoCENS)
  OdataNoCENS <- as.data.table(OdataNoCENS, key=c(ID, t))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]


  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  OData <- importData(OdataNoCENS, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")

  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR) # , gstar_TRT = ..., gstar_MONITOR = ...
  wts.DT <- getIPWeights(OData = OData, rule_name = "g.obs")
  survNP_res <- survNPMSM(wts.DT, OData)
  survNPIPW_ests <- survNP_res$IPW_estimates
  survNPIPW_ests[]

  survDirIPW_ests <- survDirectIPW(wts.DT, OData)
  survDirIPW_ests[]

  survMSM_res <- survMSM(wts.DT, OData)
  names(survMSM_res)
  survMSM_ests <- survMSM_res$St
  survMSM_ests

  res <- stremr(OdataNoCENS, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
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

  OData <- importData(OdataNoCENS, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  wts.DT <- getIPWeights(OData = OData) # , gstar_TRT = ..., gstar_MONITOR = ...)
  survNP_res <- survNPMSM(wts.DT, OData)
  survNPIPW_ests <- survNP_res$IPW_estimates

  res <- stremr(OdataNoCENS, ID = "ID", t = "t",
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
}
