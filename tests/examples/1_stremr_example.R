#-------------------------------------------------------------------
# EXAMPLE WITH CATEGORICAL CENSORING (3 levels)
#-------------------------------------------------------------------
require("data.table")
require("magrittr")
data(OdataCatCENS)
OdataDT <- as.data.table(OdataCatCENS, key=c("ID", "t"))
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

# # Estimate IPW-based hazard and survival (KM) for a rule "dhigh":
# IPW_KM_res <- stremr(OdataDT, intervened_TRT = "dhigh", intervened_MONITOR = "gstar.N",
#               ID = "ID", t = "t", covars = c("highA1c", "lastNat1"),
#               CENS = "CatC", gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
#               TRT = "TI", gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
#               MONITOR = "N", gform_MONITOR = gform_MONITOR, OUTCOME = "Y.tplus1")

# # Survival estimates by time:
# IPW_KM_res$estimates
# # Input data:
# IPW_KM_res$dataDT
