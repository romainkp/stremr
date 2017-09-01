require("stremr")
options(stremr.verbose = FALSE)
require("data.table")
data(OdataNoCENS)
Odat_DT <- as.data.table(OdataNoCENS, key=c(ID, t))

# define lagged N, first value is always 1 (always monitored at the first time point):
Odat_DT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
Odat_DT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
Odat_DT[, "continA" := rnorm(nrow(Odat_DT))] ## observed exposure
Odat_DT[, "new_continA" := continA + 0.1] ## intervened / counterfactual exposure

OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)

gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
gform_TRT = "continA ~ CVD + highA1c + N.tminus1"

test_that("Propensity score fitting with continuous exposure will fail by default", {
  # gform_MONITOR <- "N ~ 1"
  # gform_MONITOR = gform_MONITOR
  expect_error(OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT))
  OData$dat.sVar[]

  OData$modelfit.gA$predict(OData)
  options(stremr.verbose = FALSE)
})

test_that("Propensity score fitting with continuous exposure will work with Lrnr_condensier by default", {
  # options(sl3.verbose = TRUE)
  OData <- importData(Odat_DT[1:2000, ], ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData = OData,
                        gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        models_TRT = Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = FALSE))
  OData$g_preds
})

test_that("Propensity score fitting with continuous exposure will work with SuperLearner for cond. densties", {
  options(sl3.verbose = TRUE)
  options(stremr.verbose = TRUE)
  OData <- importData(Odat_DT[1:2000, ], ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
  lrn1 <- Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = FALSE)
  lrn2 <- Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = TRUE)
  sl <- Lrnr_sl$new(learners = Stack$new(lrn1, lrn2),
                    metalearner = Lrnr_solnp_density$new())
  OData <- fitPropensity(OData = OData,
                        gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        models_TRT = sl)
})

test_that("IPW works with continuous exposure and cond. dens. SuperLearner", {
  options(sl3.verbose = TRUE)
  options(stremr.verbose = TRUE)
  OData <- importData(Odat_DT[1:2000, ], ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
  lrn1 <- Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = FALSE)
  lrn2 <- Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = TRUE)
  sl <- Lrnr_sl$new(learners = Stack$new(lrn1, lrn2),
                    metalearner = Lrnr_solnp_density$new())
  OData <- fitPropensity(OData = OData,
                        gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        models_TRT = sl)
  OData$g_preds

  ## todo: now we need to somehow call sl$predict() on the counterfactual exposure as the outcome (evalute another likelihood)
  ## currently just calling ModelDeterministic$predict() which simply evaluates I(continA=new_continA), which isn't enough
  wtsDT <- getIPWeights(OData, intervened_TRT = "new_continA")
  survNPMSM <- survNPMSM(wtsDT, OData)
  survNPMSM <- survDirectIPW(wtsDT, OData)



  t.surv <- c(0:3)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)

  options(stremr.verbose = FALSE)
})


