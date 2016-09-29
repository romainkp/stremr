test.ContinModel <- function() {
  options(stremr.verbose = TRUE)
  require("data.table")
  data(OdataNoCENS)
  # data(OdatDT_10K)
  # Odat_DT <- OdatDT_10K
  Odat_DT <- as.data.table(OdataNoCENS, key=c(ID, t))

  # define lagged N, first value is always 1 (always monitored at the first time point):
  Odat_DT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  Odat_DT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
  Odat_DT[, "continA" := rnorm(nrow(Odat_DT))]

  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "continA ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  OData$dat.sVar[]

  OData$modelfit.gA$predict(OData)
  options(stremr.verbose = FALSE)
}