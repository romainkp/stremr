# test various regression / subsetting schemes and make sure it works as expected
test.DataStorageClass <- function() {
  options(stremr.verbose = TRUE)
  require("data.table")
  data(OdataNoCENS)
  OdataNoCENS <- as.data.table(OdataNoCENS, key=c(ID, t))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, "continA" := rnorm(nrow(OdataNoCENS))]
  OdataNoCENS[, "catA" := sample.int(3, size = nrow(OdataNoCENS), replace = TRUE)]

  OData <- importData(OdataNoCENS, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "continA ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  # testing various methods of DataStorageClass:
  OData$addYnode(OdataNoCENS[["Y.tplus1"]])
  checkException(OData$get.outvar(var = "blah"))

  OData$names.sVar

  # ------------------------------------------------------------------------------
  # Binning categorical variable
  # ------------------------------------------------------------------------------
  binnedCatA <- OData$binirize.sVar("catA", OData$detect.cat.sVar.levels("catA"))

  # ------------------------------------------------------------------------------
  # Binning continuous variable with different cut levels
  # ------------------------------------------------------------------------------
  intrvls1 <- OData$detect.sVar.intrvls("continA", nbins = 10, bin_bymass = FALSE, bin_bydhist = FALSE, max_nperbin = NULL)
  binnedContinA1 <- OData$binirize.sVar("continA", intervals = intrvls1, nbins = length(intrvls1)-1, bin.nms = OData$bin.nms.sVar("continA", length(intrvls1)-1))
  var.bw1 <- OData$get.sVar.bw("continA", intervals = intrvls1)
  get.sVar.bwdiff1 <- OData$get.sVar.bwdiff("continA", intervals = intrvls1)


  intrvls2 <- OData$detect.sVar.intrvls("continA", nbins = 1, bin_bymass = TRUE, bin_bydhist = FALSE, max_nperbin = 100)
  binnedContinA2 <- OData$binirize.sVar("continA", intervals = intrvls2, nbins = length(intrvls2)-1, bin.nms = OData$bin.nms.sVar("continA", length(intrvls2)-1))
  var.bw2 <- OData$get.sVar.bw("continA", intervals = intrvls2)
  get.sVar.bwdiff2 <- OData$get.sVar.bwdiff("continA", intervals = intrvls2)

  # ------------------------------------------------------------------------------
  # replace all NAs with 0:
  # ------------------------------------------------------------------------------
  OData$fixmiss_sVar()

  # ------------------------------------------------------------------------------
  # Set variable (column) types
  # ------------------------------------------------------------------------------
  OData$def.types.sVar("binary")
  OData$get.sVar.type()

  new.types <- lapply(names(OData$dat.sVar), function(x) "categor")
  names(new.types) <- names(OData$dat.sVar)
  OData$def.types.sVar(new.types)
  OData$get.sVar.type()

  OData$def.types.sVar()
  OData$get.sVar.type()

  OData$set.sVar.type("CVD", "binary")
  OData$get.sVar.type("CVD")
  checkTrue(OData$is.sVar.bin("CVD"))

  checkTrue(OData$is.sVar.cat("catA"))
  checkTrue(OData$is.sVar.cont("continA"))

  # ------------------------------------------------------------------------------
  # Setting specific column to new values
  # ------------------------------------------------------------------------------
  OData$set.sVar("continA", new.sVarVal = OdataNoCENS[["continA"]])

  # ------------------------------------------------------------------------------
  # Converting input data from long format to wide format:
  # ------------------------------------------------------------------------------
  wideDT <- OData$convert.to.wide(bslcovars = "CVD")

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

  OData$H2O.dat.sVar
  # H2O.dat.sVar = function(dat.sVar) {
  #   if (missing(dat.sVar)) {
  #     return(private$.H2O.mat.sVar)
  #   } else {
  #     assert_that(is.H2OFrame(dat.sVar))
  #     private$.H2O.mat.sVar <- dat.sVar
  #   }
  # },

  head(OData$dat.bin.sVar)
  OData$dat.bin.sVar <- OData$dat.bin.sVar

  OData$backup.savedGstarsDT
  head(OData$noNA.Ynodevals)

  checkTrue(OData$active.bin.sVar %in% "continA")
  head(OData$ord.sVar)

  # wipe out binirized .mat.sVar:
  OData$emptydat.bin.sVar
  checkTrue(is.null(OData$dat.bin.sVar))
  # wipe out input data:
  OData$emptydat.sVar
  checkTrue(is.null(OData$dat.sVar))
  options(stremr.verbose = FALSE)
}
