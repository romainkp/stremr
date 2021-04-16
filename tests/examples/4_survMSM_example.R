#-------------------------------------------------------------------
# Simulated data with informative right-censoring
#-------------------------------------------------------------------
require("data.table")
require("magrittr")
data(OdataCatCENS)
OdataDT <- as.data.table(OdataCatCENS, key=c("ID", "t"))
# Indicator that the person has never been treated in the past:
OdataDT[, "barTIm1eq0" := as.integer(c(0, cumsum(TI)[-.N]) %in% 0), by = ID]
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

#-------------------------------------------------------------------
# Import data into stremr object
#-------------------------------------------------------------------
OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"),
                    CENS = "CatC", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")

# ----------------------------------------------------------------------
# Look at the input data object
# ----------------------------------------------------------------------
print(OData)

# ----------------------------------------------------------------------
# Access the input data
# ----------------------------------------------------------------------
get_data(OData)

#-------------------------------------------------------------------
# Estimate Propensity scores
#-------------------------------------------------------------------
OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                      stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

#-------------------------------------------------------------------
# Defining weights for two dynamic regimens "dlow" and "dhigh"
#-------------------------------------------------------------------
wts.St.dlow <- getIPWeights(OData, intervened_TRT = "dlow")
wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "dhigh")

# ------------------------------------------------------------------
# Estimate survival with IPW-adjusted MSM for the hazard (logistic model)
# 1. Saturated model for time points 0 to 7
# 2. Smoothing with one hazard coefficient over time-points 8 to 11
# 3. Smoothing with one hazard coefficient over time-points 12 to 15
# ------------------------------------------------------------------
IPW_MSM_res <- survMSM(OData, wts_data = list(dlow = wts.St.dlow, dhigh = wts.St.dhigh),
                      tbreaks = c(1:8,12,16)-1,
                      est_name = "IPAW", getSEs = TRUE)
names(IPW_MSM_res)
# Survival estimates over time
IPW_MSM_res$St
# SE estimates for each time point:
IPW_MSM_res$IC.Var.S.d$dhigh$se.S
IPW_MSM_res$IC.Var.S.d$dlow$se.S
# MSM coefficient fits
IPW_MSM_res$MSM.fit

# ------------------------------------------------------------------
# Generate automatic html report with results of the analysis
# This assumes that pandoc is already installed
# For more information, go to: http://pandoc.org/installing.html
# ------------------------------------------------------------------
\dontrun{
make_report_rmd(OData, MSM = IPW_MSM_res,
  AddFUPtables = TRUE,
  RDtables = get_MSM_RDs(IPW_MSM_res, t.periods.RDs = c(12, 15), getSEs = FALSE),
  WTtables = get_wtsummary(IPW_MSM_res$wts_data,
    cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
  file.name = "sim.data.example.fup",
  title = "Custom Report Title",
  author = "Jane Doe",
  y_legend = 0.95)
}
