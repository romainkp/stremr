# ------------------------------------------------------------------------------------------------------
# SIMULATE, SAVE AND COMPRESS THE DATASET FROM THE EXAMPLE
# ------------------------------------------------------------------------------------------------------
notrun.save.example.data.02 <- function() {
  require("simcausal")
  # Nsize <- 50000
  Nsize <- 1000
  OdataCatCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356, catC=TRUE)
  OdataCatCENS[OdataCatCENS[,"t"]%in%16,"lastNat1"] <- NA
  save(OdataCatCENS, compress = TRUE, file = "./data/OdataCatCENS.rda", compression_level = 9)
  require("tools")
  resaveRdaFiles("./data/OdataCatCENS.rda", compress = "bzip2")

  # --------------------------------
  # save data as csv
  # --------------------------------
  # data.table::fwrite(OdataCatCENS, "./OdataCatCENS.csv", turbo = TRUE, verbose = TRUE, na = "NA_h2o")
}

test.model.fits.categorCENSOR <- function() {
  require("data.table")
  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  data(OdataCatCENS)
  OdataCatCENS <- as.data.table(OdataCatCENS, key=c(ID, t))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataCatCENS[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataCatCENS[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]

  # head(OdataCatCENS)
  # nrow(OdataCatCENS)
  # table(OdataCatCENS[,"Y.tplus1"])
  # table(OdataCatCENS[,"C"])
  # unique(OdataCatCENS[,"C"])
  # unique(OdataCatCENS[,"CatC"])
  # [1]  0  2  1 NA

  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  res <- stremr(OdataCatCENS, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
          # noCENScat = 0L)
  res$IPW_estimates
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

  res <- stremr(OdataCatCENS, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
          gform_TRT = gform_TRT,
          gform_MONITOR = gform_MONITOR)
          # noCENScat = 0L)
  res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 3:
  # --------------------------------
  # options(stremr.verbose = TRUE)
  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  res <- stremr(OdataCatCENS, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
          gform_TRT = gform_TRT,
          gform_MONITOR = gform_MONITOR,
          noCENScat = 1L)
    # )
  res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 4:
  # --------------------------------

  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify_CENS <- list(C = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(TI=c("t == 0L", "t > 0L & TI.tminus1 == 0L", "t > 0L & TI.tminus1 == 1L"))
  # **** really want to define it like this ****
  # gform_TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
  #               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))
  gform_MONITOR <- "N ~ 1"

  # system.time(
  res <- stremr(OdataCatCENS, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
          gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
          gform_MONITOR = gform_MONITOR,
          noCENScat = 0L)
    # )
  res$IPW_estimates
  # res$dataDT

}
