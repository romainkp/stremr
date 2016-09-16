test.model.fits.categorCENSOR <- function() {
  options(stremr.verbose = TRUE)
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

test.model.fits.categorCENSOR2 <- function() {
  options(stremr.verbose = TRUE)
  #-------------------------------------------------------------------
  # EXAMPLE WITH CATEGORICAL CENSORING (3 levels)
  #-------------------------------------------------------------------
  require("data.table")
  require("magrittr")
  data(OdataCatCENS)
  OdataDT <- as.data.table(OdataCatCENS, key=c(ID, t))
  # Indicator that the person has never been treated in the past:
  OdataDT[, "barTIm1eq0" := as.integer(c(0, cumsum(TI)[-.N]) %in% 0), by = ID]
  # Define lagged N, first value is always 1 (always monitored at the first time point):
  OdataDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]

  #-------------------------------------------------------------------
  # Regressions for modeling the exposure (TRT)
  #-------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  # Fit a separate model for TRT (stratify) for each of the following subsets:
  stratify_TRT <- list(
    TI=c(
         # MODEL TI AT t=0
         "t == 0L",
         # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",
         # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",
         # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
         "(t > 0L) & (barTIm1eq0 == 1L)"
        ))

  #-------------------------------------------------------------------
  # Regressions for modeling the categorical censoring (CENS)
  #-------------------------------------------------------------------
  gform_CENS <- c("CatC ~ highA1c")
  # stratify by time-points (separate model for all t<16 and t=16)
  stratify_CENS <- list(CatC=c("t < 16", "t == 16"))

  #-------------------------------------------------------------------
  # Regressions for modeling the monitoring regimen (MONITOR)
  #-------------------------------------------------------------------
  # Intercept only model, pooling across all time points t
  gform_MONITOR <- "N ~ 1"

  #-------------------------------------------------------------------
  # Define the counterfactual monitoring regimen of interest
  #-------------------------------------------------------------------
  # probability of being monitored at each t is 0.1
  OdataDT[, "gstar.N" := 0.1]

  # Define two dynamic rules: dlow & dhigh
  OdataDT <- defineIntervedTRT(OdataDT, theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
                            CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1",
                            new.TRT.names = c("dlow", "dhigh"), return.allcolumns = TRUE)

  # Estimate IPW-based hazard and survival (KM) for a rule "dhigh":
  res <- stremr(OdataDT, intervened_TRT = "dhigh", intervened_MONITOR = "gstar.N",
                ID = "ID", t = "t", covars = c("highA1c", "lastNat1"),
                CENS = "CatC", gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
                TRT = "TI", gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
                MONITOR = "N", gform_MONITOR = gform_MONITOR, OUTCOME = "Y.tplus1")

  res$IPW_estimates
  res$dataDT
}