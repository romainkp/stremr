# test various regression / subsetting schemes and make sure it works as expected
test.DataStorageClass <- function() {
  options(stremr.verbose = TRUE)
  require("data.table")
  data(OdataNoCENS)
  OdataNoCENS <- as.data.table(OdataNoCENS, key=c("ID", "t"))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, "continA" := rnorm(nrow(OdataNoCENS))]
  OdataNoCENS[, "catA" := sample.int(3, size = nrow(OdataNoCENS), replace = TRUE)]

  OData <- importData(OdataNoCENS, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "continA ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  ## testing the printing functionality and testing the functionality to obtain the data
  OData
  OData$dat.sVar
  print(OData)

  intputDT <- get_data(OData)
  intputDT[]

  # testing various methods of DataStorageClass:
  OData$addYnode(OdataNoCENS[["Y.tplus1"]])
  checkException(OData$get.outvar(var = "blah"))

  OData$names.sVar

  # ------------------------------------------------------------------------------
  # Active bindings
  # ------------------------------------------------------------------------------
  OData$min.t
  OData$max.t
  OData$nobs
  OData$nuniqueIDs
  OData$nuniquets
  OData$names.sVar
  OData$ncols.sVar
  OData$dat.sVar


  OData$backup.savedGstarsDT
  head(OData$noNA.Ynodevals)
  options(stremr.verbose = FALSE)
}
