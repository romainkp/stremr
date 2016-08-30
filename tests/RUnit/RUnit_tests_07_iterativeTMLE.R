`%+%` <- function(a, b) paste0(a, b)
require("data.table")
require("stremr")
data(OdatDT_10K)
Odat_DT <- OdatDT_10K
# ---------------------------------------------------------------------------
# Define some summaries (lags C[t-1], A[t-1], N[t-1])
# ---------------------------------------------------------------------------
ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
lagnodes <- c("C", "TI", "N")
newVarnames <- lagnodes %+% ".tminus1"
Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
# indicator that the person has never been on treatment up to current t
Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

# ----------------------------------------------------------------
# IMPORT DATA
# ----------------------------------------------------------------
options(stremr.verbose = TRUE)
set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")
OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
# ------------------------------------------------------------------
# Fit propensity scores for Treatment, Censoring & Monitoring
# ------------------------------------------------------------------
gform_TRT <- c("TI ~ CVD + highA1c")
stratify_TRT <- list(
  TI=c("t == 0L",                                            # MODEL TI AT t=0
       "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
       "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
       "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
      ))

gform_CENS <- c("C ~ highA1c + t")
gform_MONITOR <- "N ~ 1"

OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

# ---------------------------------------------------------------------------------------------------------
# Iterative TMLE
# ---------------------------------------------------------------------------------------------------------
t.surv <- c(10)
Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
params = list(fit.package = "speedglm", fit.algorithm = "glm")
# params = list(fit.package = "h2o", fit.algorithm = "RF", ntrees = 100,
#               learn_rate = 0.05, sample_rate = 0.8,
#               col_sample_rate = 0.8, balance_classes = TRUE)

iterTMLE_est1a <- fitIterTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
iterTMLE_est1a[]
iterTMLE_est1b <- fitIterTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
iterTMLE_est1b[]

iterTMLE_est2a <- fitIterTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
iterTMLE_est2a[]
iterTMLE_est2b <- fitIterTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
iterTMLE_est2b[]

