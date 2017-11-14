context("CV TMLE")

reqxgb <- requireNamespace("xgboost", quietly = TRUE)
reqh2o <- requireNamespace("h2o", quietly = TRUE)
if (!reqxgb) {
# if (!reqxgb || !reqh2o) {
  run_test <- FALSE
} else {
  run_test <- TRUE
  library("xgboost")
  # library("h2o")
}

library("stremr")
library("magrittr")
library("data.table")
library("testthat")
library("sl3")

data.table::setDTthreads(1)
options(stremr.verbose = FALSE)
options(gridisl.verbose = FALSE)
options(sl3.verbose = FALSE)
options(condensier.verbose = FALSE)
# options(stremr.verbose = TRUE)
# options(gridisl.verbose = TRUE)
# options(sl3.verbose = TRUE)
# options(condensier.verbose = TRUE)


data(OdatDT_10K)
Odat_DT <- OdatDT_10K
Odat_DT <- Odat_DT[ID %in% (1:100), ] ## select only the first 100 IDs (100 subjects)
setkeyv(Odat_DT, cols = c("ID", "t"))

## ---------------------------------------------------------------------------
## Define some summaries (lags C[t-1], A[t-1], N[t-1])
## ---------------------------------------------------------------------------
ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
lagnodes <- c("C", "TI", "N")
newVarnames <- paste0(lagnodes, ".tminus1")
Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
## indicator that the person has never been on treatment up to current t
Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

## ------------------------------------------------------------------
## Propensity score models for Treatment, Censoring & Monitoring
## ------------------------------------------------------------------
gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
stratify_TRT <- list(
  TI=c("t == 0L",                                            # MODEL TI AT t=0
       "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
       "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
       "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
      ))

gform_CENS <- c("C ~ highA1c + t")
gform_MONITOR <- "N ~ 1"

## ----------------------------------------------------------------
## IMPORT DATA
## ----------------------------------------------------------------
OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
OData$dat.sVar[]
OData$fold_column <- NULL
OData$nfolds <- NULL

## ----------------------------------------------------------------
## FIT PROPENSITY SCORES WITH xgboost gbm and V fold CV
## ----------------------------------------------------------------
## xgboost gbm
params_g <- gridisl::defModel(estimator = "xgboost__gbm", family = "quasibinomial", rounds = 5, nthread = 1)
OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                        stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                        models_CENS = params_g, models_TRT = params_g, models_MONITOR = params_g,
                        fit_method = "cv", fold_column = "fold_ID")
## h2o gbm
# OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
#                        stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
#                        estimator = "h2o__gbm", distribution = "bernoulli",
#                        models_MONITOR = defModel(estimator = "speedglm__glm", family = "quasibinomial"),
#                       fit_method = "cv", fold_column = "fold_ID", ntrees = 5
#                       )

wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
surv_dlow <- survNPMSM(wts.St.dlow, OData) %$% estimates
wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
surv_dhigh <- survNPMSM(wts.St.dhigh, OData) %$% estimates

## weights based on holdout propensity score predictions
wts.St.dlow_hold <- getIPWeights(OData, intervened_TRT = "gTI.dlow", holdout = TRUE)
surv_dlow_hold <- survNPMSM(wts.St.dlow_hold, OData) %$% estimates
wts.St.dhigh_hold <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", holdout = TRUE)
surv_dhigh_hold <- survNPMSM(wts.St.dhigh_hold, OData) %$% estimates

test_that("regular and holdout CV IPW estimates are as expected", {
  if (run_test) {
    expect_equal(
      surv_dlow[["St.NPMSM"]],
      c(0.9857093, 0.9857093, 0.9857093, 0.9857093, 0.9857093, 0.9714186, 0.9714186, 0.8589230, 0.8589230, 0.8589230, 0.8589230, 0.8589230, 0.8446323, 0.8446323, 0.8303416, 0.8303416, NaN),
      tolerance = .0001
    )
    expect_equal(
      surv_dlow_hold[["St.NPMSM"]],
      c(0.9874289, 0.9874289, 0.9874289, 0.9874289, 0.9874289, 0.9694629, 0.9694629, 0.8809416, 0.8809416, 0.8809416, 0.8809416, 0.8809416, 0.8629710, 0.8629710, 0.8449992, 0.8449992, NaN),
      tolerance = .0001
    )
    expect_equal(
      surv_dhigh[["St.NPMSM"]],
      c(0.9871938, 0.9774610, 0.9774610, 0.9310049, 0.8697201, 0.8358002, 0.7917383, 0.7330401, 0.7330401, 0.6906911, 0.6906911, 0.6410277, 0.6410277, 0.6410277, 0.6410277, 0.6410277, NaN),
      tolerance = .0001
    )
    expect_equal(
      surv_dhigh_hold[["St.NPMSM"]],
      c(0.9878905, 0.9782232, 0.9782232, 0.9102207, 0.8502852, 0.8183632, 0.7766260, 0.7178761, 0.7178761, 0.6764124, 0.6764124, 0.6261653, 0.6261653, 0.6261653, 0.6261653, 0.6261653, NaN),
      tolerance = .0001
    )
  }
})

# ---------------------------------------------------------------------------------------------------------
# CV TMLE w/ xgboost gbm and cross-validation selection of Q
# ---------------------------------------------------------------------------------------------------------
params <- defModel(estimator = "xgboost__gbm",
                    family = "quasibinomial",
                    nthread = 1,
                    nrounds = 5)

t.surv <- c(0:10)
Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

CV_tmle_est <- fit_CVTMLE(OData, tvals = t.surv,
                          intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                          stratifyQ_by_rule = FALSE,
                          fit_method = "cv",
                          fold_column = "fold_ID",
                          parallel = FALSE)

test_that("CV_TMLE point estimates and SEs are as expected", {
  if (run_test) {
    expect_equal(
      CV_tmle_est[["estimates"]][["St.TMLE"]],
      c(0.9884258, 0.9791671, 0.8619628, 0.9126059, 0.8538316, 0.8197174, 0.7769200, 0.7192738, 0.6317145, 0.6788607, 0.5962425),
      tolerance = .0001
    )
    expect_equal(
      CV_tmle_est[["estimates"]][["SE.TMLE"]],
      c(0.01405839, 0.01805242, 0.02968440, 0.06479710, 0.07355126, 0.07649550, 0.08180215, 0.08693926, 0.08957773, 0.09119548, 0.09393821),
      tolerance = .0001
    )
  }
})

tmle_est <- fit_TMLE(OData, tvals = t.surv,
                    intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                    stratifyQ_by_rule = FALSE,
                    fit_method = "cv",
                    fold_column = "fold_ID",
                    parallel = FALSE)

test_that("TMLE estimates with CV are as expected", {
  if (run_test) {
    expect_equal(
      tmle_est[["estimates"]][["St.TMLE"]],
      c(0.9872006, 0.9774787, 0.8629965, 0.9310966, 0.8704412, 0.8330573, 0.7880491, 0.7275443, 0.6391764, 0.6860308, 0.6067204),
      tolerance = .0001
    )
    expect_equal(
      tmle_est[["estimates"]][["SE.TMLE"]],
      c(0.01424316, 0.01800686, 0.02809040, 0.04373658, 0.05765227, 0.06225333, 0.06978813, 0.07720831, 0.08191416, 0.08218477, 0.08668708),
      tolerance = .0001
    )
  }
})

