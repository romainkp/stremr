#-------------------------------------------------------------------
# EXAMPLE BASED ON SIMULATED DATA
#-------------------------------------------------------------------
library("data.table")
library("magrittr")
data(OdataCatCENS)
OdataDT <- as.data.table(OdataCatCENS, key=c(ID, t))
# Indicator that the person has never been treated in the past:
OdataDT[, "barTIm1eq0" := as.integer(c(0, cumsum(TI)[-.N]) %in% 0), by = ID]

#-------------------------------------------------------------------
# Regressions for modeling the exposure (TRT)
#-------------------------------------------------------------------
gform.TRT <- "TI ~ CVD + highA1c + N.tminus1"
# Fit a separate model for TRT (stratify) for each of the following subsets:
stratify.TRT <- list(
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
gform.CENS <- c("CatC ~ highA1c")
# stratify by time-points (separate model for all t<16 and t=16)
stratify.CENS <- list(CatC=c("t < 16", "t == 16"))

#-------------------------------------------------------------------
# Regressions for modeling the monitoring regimen (MONITOR)
#-------------------------------------------------------------------
# Intercept only model, pooling across all time points t
gform.MONITOR <- "N ~ 1"

#-------------------------------------------------------------------
# Define the counterfactual monitoring regimen of interest
#-------------------------------------------------------------------
p <- 0.1 # probability of being monitored at each t is 0.1
OdataDT[, "gstar.N" := ifelse(N == 1L, eval(p), 1-eval(p))]

# Define rule followers/non-followers for two rules: dlow & dhigh
res <- follow.rule.d.DT(OdataDT,
        theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
        CENS = "CatC", TRT = "TI", MONITOR = "N",
        rule.names = c("dlow", "dhigh")) %>%
# Merge rule definitions into main dataset:
  merge(OdataDT, ., by=c("ID", "t")) %>%
# Estimate hazard and survival for a rule "dhigh":
