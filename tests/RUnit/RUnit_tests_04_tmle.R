`%+%` <- function(a, b) paste0(a, b)
# --------------------------------
# INSTALL CORRECT VERSIONS of data.table and stremr from github:
# --------------------------------
# devtools::install_github('Rdatatable/data.table')
require("data.table")
# devtools::install_github('osofr/stremr', ref = "tmle", build_vignettes = FALSE)
require("stremr")

# --------------------------------
# Test data set included in stremr:
# --------------------------------
data(O.data.simstudy.g05)
O.data <- O.data.simstudy.g05
head(O.data)

ID <- "ID"; t <- "t"; TRT <- "TI"; CENS <- "C"; MONITOR <- "N"; outcome <- "Y"; I <- "highA1c";

# --------------------------------
# Define counterfactual treatment assignment under two rules (dlow & dhigh)
# --------------------------------
O.dataDT <- data.table(O.data, key = c(ID, t))
# Counterfactual TRT assignment for rule dlow (equivalent to always treated):
rule_name1 <- "dlow"
O.dataDT[,"TI.gstar." %+% rule_name1 := 1L]
# Counterfactual TRT assignment for dynamic rule dhigh -> start TRT only when I=1 (highA1c = 1)
rule_name2 <- "dhigh"
O.dataDT_TIdhigh <- stremr::defineIntervedTRT(O.dataDT, theta = 1, ID = ID,
                                          t = t, I = I, CENS = CENS, TRT = TRT,
                                          MONITOR = MONITOR,
                                          tsinceNis1 = "lastNat1",
                                          new.TRT.names = "TI.gstar." %+% rule_name2)
O.dataDT <- merge(O.dataDT, O.dataDT_TIdhigh, by=c(ID, t))

# ---------------------------------------------------------------------------
# DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
# Might expand this in the future to allow defining arbitrary summaries
# ---------------------------------------------------------------------------
lagnodes <- c("C", "TI", "N")
newVarnames <- lagnodes %+% ".tminus1"
O.dataDT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]

# indicator that the person has never been on treatment up to current t
TIcovarname <- "barTIm1eq0"
O.dataDT[, (TIcovarname) := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]

# -------------------------------------------------------------------------------------------
# Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
# -------------------------------------------------------------------------------------------
OUTCOME <- "Y"
shifted.OUTCOME <- OUTCOME%+%".tplus1"
O.dataDT[, (shifted.OUTCOME) := shift(get(OUTCOME), n = 1L, type = "lead"), by = ID]
O.dataDT <- O.dataDT[!get(OUTCOME)%in%1,]

O.dataDT[1:100, ]

# --------------------------------
# Define global options for stremr (which R packages to use for model fitting)
# --------------------------------
# options(stremr.verbose = FALSE)
options(stremr.verbose = TRUE)
stremr_options(fit.package = "speedglm", fit.algorithm = "GLM")

# import data into stremr object:
OData <- importData(O.dataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = shifted.OUTCOME)

# --------------------------------
# Fitting the propensity scores for observed variables (A,C,N)
# --------------------------------
# + N.tminus1
gform_TRT <- "TI ~ CVD + highA1c"
stratify_TRT <- list(
  TI=c("t == 0L",                                            # MODEL TI AT t=0
       "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
       "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
       "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
      ))
gform_CENS <- c("C ~ highA1c")
stratify_CENS <- list(C=c("t < 16", "t == 16"))
gform_MONITOR <- "N ~ 1"

OData <- fitPropensity(OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT,
                              stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

# get IPW-adjusted and KM survival (with hazards over time)
wts.St.dlow <- getIPWeights(OData, intervened_TRT = "TI.gstar.dlow")
St.dlow <- survNPMSM(wts.St.dlow, OData)
St.dlow

wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "TI.gstar.dhigh")
St.dhigh <- survNPMSM(wts.St.dhigh, OData)
St.dhigh

# ---------------------------------------------------------------------------------------------------------
# GCOMP AND TMLE w/ GLMs
# ---------------------------------------------------------------------------------------------------------
t.surv <- c(1,2,3,4,5,6,7,8,9,10)
# t.surv <- c(1,2,3)

Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

# stratified modeling by rule followers only:
gcomp_est1 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
tmle_est1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
gcomp_est1; tmle_est1

# pooling all observations (no stratification):
gcomp_est2 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
tmle_est2 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
gcomp_est2; tmle_est2

# stratified modeling by rule followers only:
gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
gcomp_est3; tmle_est3

# pooling all observations (no stratification):
gcomp_est4 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
gcomp_est4; tmle_est4

# ------------------------------------------------------------------------
# RUN PARALLEL seq-GCOMP & TMLE over t.surv (MUCH FASTER)
# ------------------------------------------------------------------------
require("doParallel")
registerDoParallel(cores = 4)
data.table::setthreads(1)

gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
gcomp_est; tmle_est

gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
gcomp_est; tmle_est

# ---------------------------------------------------------------------------------------------------------
# GCOMP AND TMLE w/ h2o random forest
# ---------------------------------------------------------------------------------------------------------
require("h2o")
h2o::h2o.init(nthreads = 4)
# h2o::h2o.init()
params = list(fit.package = "h2o", fit.algorithm = "RF", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)

t.surv <- c(1,2,3,4,5,6,7,8,9,10)
Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
gcomp_est; tmle_est

gcomp_fit <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.gstar.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
gcomp_est; tmle_est

# ------------------------------------------------------------------------
# TEST FOR ERROR WITH > 1 REGRESSION AND > 1 STATA
# ------------------------------------------------------------------------
gform_CENS_test <- c("C1 ~ highA1c", "C2 ~ highA1c")
stratify_CENS_test <- list(C1=c("t < 16", "t == 16"), C2=c("t < 16", "t == 16"))
O.dataDT_test <- O.dataDT
O.dataDT_test[, "C1" := C]
O.dataDT_test[, "C2" := C]
OData <- importData(O.dataDT_test, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = c("C1","C2"), TRT = "TI", MONITOR = "N", OUTCOME = shifted.OUTCOME)
OData <- fitPropensity(OData, gform_CENS = gform_CENS_test, stratify_CENS = stratify_CENS_test, gform_TRT = gform_TRT,
                              stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                              params_CENS = params_CENS, params_TRT = params_TRT, params_MONITOR = params_MONITOR)
