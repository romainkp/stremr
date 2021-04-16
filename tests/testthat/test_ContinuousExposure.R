context("continuous exposure")

require("stremr")
require("sl3")
require("data.table")

options(stremr.verbose = FALSE)
options(gridisl.verbose = FALSE)
options(sl3.verbose = FALSE)
options(condensier.verbose = FALSE)
# options(stremr.verbose = TRUE)
# options(gridisl.verbose = TRUE)
# options(sl3.verbose = TRUE)
# options(condensier.verbose = TRUE)

data(OdataNoCENS)
Odat_DT <- as.data.table(OdataNoCENS, key=c("ID", "t"))

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
  # expect_error()
  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT)
  OData$dat.sVar[]

  OData$modelfit.gA$predict(OData)
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

test_that("IPW and TMLE run with continuous exposure and cond. dens. SuperLearner", {
  OData <- importData(Odat_DT[1:2000, ], ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "continA", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
  lrn1 <- Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = FALSE)
  lrn2 <- Lrnr_condensier$new(nbins = 5, bin_method = "equal.mass", pool = TRUE)
  sl <- Lrnr_sl$new(learners = Stack$new(lrn1, lrn2),
                    metalearner = Lrnr_solnp_density$new())

  stratify_TRT <- list(continA=c("t == 0L", "t > 0L"))

  OData <- fitPropensity(OData = OData,
                         gform_CENS = gform_CENS,
                         gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT,
                         models_TRT = sl)
  # OData$g_preds

  ## todo: now we need to somehow call sl$predict() on the counterfactual exposure as the outcome (evalute another likelihood)
  ## currently just calling ModelDeterministic$predict() which simply evaluates I(continA=new_continA), which isn't enough
  wtsDT <- getIPWeights(OData, intervened_TRT = "new_continA", intervened_type_TRT = "shift")
  survNPMSM <- survNPMSM(wtsDT, OData)
  survMSM <- survMSM(wtsDT, OData)
  IPW <- directIPW(wtsDT, OData)

  ## gstar.A should be 1
  wtsDT_MSM <- getIPWeights(OData, intervened_TRT = "new_continA", intervened_type_TRT = "MSM")

  t.surv <- c(0:3)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "new_continA", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "new_continA",
                        Qforms = Qforms,
                        stratifyQ_by_rule = FALSE,
                        intervened_type_TRT = "shift")

  options(stremr.verbose = FALSE)
})


