run_test <- TRUE
library("stremr")
library("sl3")
library("magrittr")
library("data.table")
library("testthat")

data.table::setDTthreads(1)
options(stremr.verbose = TRUE)
options(gridisl.verbose = TRUE)
options(sl3.verbose = TRUE)

## --------------------------------------------------------------------------------------------
## Perform estimation with LTMLE for stochastic intervention g^*
## Using simulated observed data
## --------------------------------------------------------------------------------------------
## load 2 time-point data example
data(dt_2t)
true_EYgstar <- attributes(dt_2t)[["true_EYgstar"]]
true_EYA0 <- attributes(dt_2t)[["true_EYA0"]]
true_EYA1 <- attributes(dt_2t)[["true_EYA1"]]

## define the counterfactual A^*
dt_2t[, ("Astoch") := 0.3][, ("A0") := 0][, ("A1") := 1]
## add some lags
lagnodes <- c("A", "L1", "L2", "L3")
newVarnames <- paste0(lagnodes, ".tminus1")
dt_2t[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
OData <- importData(dt_2t, ID = "ID", t = "t", covars = c("L1", "L2", "L3", "A.tminus1", "L1.tminus1", "L2.tminus1", "L3.tminus1"), TRT = "A", OUTCOME = "Y")
OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID")
## -------------------------------------------------
## CORRECT g
## Check that IPW is unbiased
## -------------------------------------------------
gform_TRT = "A ~ L1 + L2 + L3 + A.tminus1"
stratify_TRT <- list(A=c("t == 0", "t == 1"))

gridisl_cv_t <- system.time({
  models_TRT <- defModel(estimator = "xgboost__gbm", family = "quasibinomial", nrounds = 5)
  OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = models_TRT, fit_method = "cv")
})

sl3_cv_t <- system.time({
  lrn_xgb <- Lrnr_xgboost$new(nrounds = 5, outcome_type = "binomial")
  lrn_glm <- Lrnr_glm_fast$new(outcome_type = "binomial")
  sl <- Lrnr_sl$new(learners = Stack$new(lrn_xgb, lrn_glm), metalearner = Lrnr_solnp$new())
  OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = sl)
})

# models_TRT <- defModel(estimator = "speedglm__glm", family = "quasibinomial")
# OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = models_TRT)

# sl3_cv_h2o_t <- system.time({
#   lrn_h2o <- Lrnr_h2o_grid$new(algorithm = "gbm", ntrees = 20, outcome_type = "binomial", family = "bernoulli")
#   sl <- Lrnr_sl$new(learners = Stack$new(lrn_h2o), metalearner = Lrnr_nnls$new())
#   OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = sl)
# })

# gridisl_cv_h2o_t <- system.time({
#   # models_TRT <- defModel(estimator = "speedglm__glm")
#   models_TRT <- defModel(estimator = "h2o__gbm", ntrees = 20, distribution = "bernoulli")
#   OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = models_TRT, fit_method = "cv")
# })

# print("sl3_t"); print(sl3_t)
 #   user  system elapsed
 # 30.464   1.346  31.749
# print("gridisl_t"); print(gridisl_t)
# [1] "gridisl_t"
#    user  system elapsed
#   2.639   0.906   3.549

IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYgstar <- 1-IPW.St[time == 1, ][, St.directIPW])
cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
## [1] 0.5605612
# IPW bias g^* Astoch:  -0.0006032359

IPW.St <- getIPWeights(OData, intervened_TRT = "A0") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYA0 <- 1-IPW.St[time == 1, ][, St.directIPW])
cat("\nIPW bias A0: ", true_EYA0-IPW_EYA0, "\n")
## [1] 0.4987714
## IPW bias A0:  0.001484633

IPW.St <- getIPWeights(OData, intervened_TRT = "A1") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYA1 <- 1-IPW.St[time == 1, ][, St.directIPW])
cat("\nIPW bias A1: ", true_EYA1-IPW_EYA1, "\n")
## [1] 0.7664964
## IPW bias A1:  -0.0007333567

test_that("Stochastic g^*: IPW unbiased for correct g", {
  if (run_test) {
    expect_equal(true_EYgstar - IPW_EYgstar, -0.000604054)
    expect_equal(true_EYA0 - IPW_EYA0, 0.00148446)
    expect_equal(true_EYA1 - IPW_EYA1, -0.000733655)
  }
})

## -------------------------------------------------
## CORRECT g and CORRECT Q
## 1. sequential GCOMP
## 2. TMLE
## -------------------------------------------------
Qforms <- rep.int("Qkplus1 ~ L1 + L2 + L3 + A + A.tminus1 + L1.tminus1 + L2.tminus1 + L3.tminus1", 2)

## Define interaction terms as a separate learner:
lnr_inter <- Lrnr_define_interactions$new(interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
## Define Super-Learner candidate libraries:
lrn_xgb <- Lrnr_xgboost$new(outcome_type = "quasibinomial", nrounds = 20)
lrn_glm <- Lrnr_glm_fast$new(outcome_type="quasibinomial")
lrn_glm2 <- Lrnr_glm_fast$new(outcome_type="quasibinomial", covariates = c("A", "L1", "L2", "L3"))
lrn_glmnet <- Lrnr_glmnet$new(outcome_type = "quasibinomial", nlambda = 5)
lrn_stack <- Stack$new(lrn_xgb, lrn_glm, lrn_glm2, lrn_glmnet) ## Stack candidates:
lrn_pipe <- Pipeline$new(lnr_inter, lrn_stack) ## Pipeline interaction terms into stack of learners
lrn_sl <- Lrnr_sl$new(learners = lrn_stack, metalearner = Lrnr_solnp$new()) ## Define the super-learner on the entire thing with metalearner
# Qmodels <- Lrnr_sl$new(learners = list(lrn_pipe), metalearner = Lrnr_solnp$new()) ## Define the super-learner on the entire thing with metalearner
Qmodels <- lrn_sl
# Qmodels <- Pipeline$new(lnr_inter, lrn_sl)

# Qmodels <- gridisl::defModel(estimator = "xgboost__gbm", nrounds = 20, interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
gcomp_est <- fit_GCOMP(OData, tvals = 1, intervened_TRT = "Astoch", Qforms = Qforms, models = Qmodels)
(GCOMP_EYgstar <- 1 - gcomp_est$estimates[, St.GCOMP])
cat("\nGCOMP bias g^* Astoch: ", true_EYgstar-GCOMP_EYgstar, "\n")
## (glm) [1] 0.5586502 -- unbiased
## GCOMP bias g^* Astoch:  0.00130783

tmle_est <- fit_TMLE(OData, tvals = 1, intervened_TRT = "Astoch", Qforms = Qforms, models = Qmodels)
(TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
## [1] 0.5605454 -- unbiased
## TMLE bias g^* Astoch:  -0.0005873876

test_that("Stochastic g^*: GCOMP and TMLE are unbiased for correct g & correct Q", {
  if (run_test) {
    expect_equal(true_EYgstar - GCOMP_EYgstar, 0.001988167)
    expect_equal(true_EYgstar - TMLE_EYgstar,  -0.0005547971)
    expect_equal(tmle_est$estimates[["SE.TMLE"]], 0.0009420197)
  }
})