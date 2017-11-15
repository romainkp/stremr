run_test <- TRUE

## --------------------------------------------------------------------------------------------
## Simulate some data w/ 2 time-points
## Evaluate the truth under stochastic intervention g^* (also under two static interventions)
## --------------------------------------------------------------------------------------------
sim_stochastic_test_data <- function() {
  library("data.table")
  library("simcausal")

  set.seed(123456)

  D <- DAG.empty()
  D <- D +
    node("L1", t = 0, distr = "rbern", prob = 0.5) +
    node("L2", t = 0, distr = "rbern", prob = 0.5) +
    node("L3", t = 0, distr = "rbern", prob = 0.5) +
    node("A",  t = 0, distr = "rbern", prob = plogis(0.5*L1[0] + 0.5*L2[0] + 0.5*L3[0])) +
    node("Y",  t = 0, distr = "rbern", prob = 0) +
    node("L1", t = 1, distr = "rbern", prob = plogis(-0.5 + 0.5*L1[0] + 1.5*A[0] + 0.5*L1[0])) +
    node("L2", t = 1, distr = "rbern", prob = plogis(-0.5 + 0.5*L1[0] + 1.5*A[0] + 0.5*L2[0])) +
    node("L3", t = 1, distr = "rbern", prob = plogis(-0.5 + 0.5*L1[0] + 1.5*A[0] + 0.5*L3[0])) +
    node("A",  t = 1, distr = "rbern", prob = plogis(+1.5*A[0] - 0.5*L1[1] - 0.5*L2[1] - 0.5*L3[1])) +
    node("Y",  t = 1, distr = "rbern", prob = plogis(0.5*L1[1]*A[1] + 0.5*L2[1]*A[1] + 0.5*L3[1]*A[1]))

  Dset <- set.DAG(D)
  dt_2t <- data.table(sim(Dset, n = 500000, wide = FALSE))
  dt_2t[, lapply(.SD, mean), by = t]

  Dset <- Dset +
    action("A0", nodes = node("A", distr = "rbern", prob = 0, t = c(0:1))) +
    action("A1", nodes = node("A", distr = "rbern", prob = 1, t = c(0:1))) +
    action("Astoch", nodes = node("A", distr = "rbern", prob = 0.3, t=0:1))

  Pstar.Astoch <- data.table(sim(Dset, actions = "Astoch", n = 1000000, wide = FALSE)[[1]])
  Pstar.A0 <- data.table(sim(Dset, actions = "A0", n = 1000000, wide = FALSE)[[1]])
  Pstar.A1 <- data.table(sim(Dset, actions = "A1", n = 1000000, wide = FALSE)[[1]])

  Pstar.Astoch[, lapply(.SD, mean), by = t]
  Pstar.A0[, lapply(.SD, mean), by = t]
  Pstar.A1[, lapply(.SD, mean), by = t]

  (true_EYgstar <- mean(Pstar.Astoch[t==1, Y]))
  ## [1] 0.559198
  (true_EYA0 <- mean(Pstar.A0[t==1, Y]))
  ## [1] 0.500297
  (true_EYA1 <- mean(Pstar.A1[t==1, Y]))
  ## [1] 0.765463

  ## **** UNCOMMENT TO SAVE THIS DATASET ****
  ## WILL SAVE TO CURRENT WORKING DIRECTORY
  # notrun.save.example.data.05 <- function(dt_2t) {
  #   require("tools")
  #   save(dt_2t, compress = TRUE, file = "./dt_2t.rda", compression_level = 9)
  #   resaveRdaFiles("./dt_2t.rda", compress = "bzip2")
  # }

  # attributes(dt_2t)[["true_EYgstar"]] <- true_EYgstar
  # attributes(dt_2t)[["true_EYA0"]] <- true_EYA0
  # attributes(dt_2t)[["true_EYA1"]] <- true_EYA1
  # notrun.save.example.data.05(dt_2t)
}

library("stremr")
library("magrittr")
library("data.table")
library("testthat")

data.table::setDTthreads(1)
options(stremr.verbose = FALSE)
options(gridisl.verbose = FALSE)
options(sl3.verbose = FALSE)
options(condensier.verbose = FALSE)
# options(stremr.verbose = TRUE)
# options(gridisl.verbose = TRUE)
# options(sl3.verbose = TRUE)
# options(condensier.verbose = TRUE)

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

# models_TRT <- defModel(estimator = "speedglm__glm")
OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT)

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
    expect_equal(true_EYgstar - IPW_EYgstar, -0.0006032359)
    expect_equal(true_EYA0 - IPW_EYA0, 0.001484633)
    expect_equal(true_EYA1 - IPW_EYA1, -0.0007333567)
  }
})

## -------------------------------------------------
## CORRECT g and WRONG Q
## 1. GCOMP for stochastic intervention
## 2. TMLE
## -------------------------------------------------
Qforms <- rep.int("Qkplus1 ~ L1 + L2 + L2.tminus1 + L3.tminus1", 2)

# params <- gridisl::defModel(estimator = "speedglm__glm")
gcomp_est <- fit_GCOMP(OData, tvals = 1, intervened_TRT = "Astoch", Qforms = Qforms)
(GCOMP_EYgstar <- 1 - gcomp_est$estimates[, St.GCOMP])
cat("\nGCOMP bias g^* Astoch: ", true_EYgstar-GCOMP_EYgstar, "\n")
## [1] 0.614744 -- very biased
## GCOMP bias g^* Astoch:  -0.054786

tmle_est <- fit_TMLE(OData, tvals = 1, intervened_TRT = "Astoch", Qforms = Qforms)
(TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
## 0.5605559 -- unbiased, same as IPW
## TMLE bias g^* Astoch:  -0.0005979065

test_that("Stochastic g^*: GCOMP is biased and TMLE is unbiased for correct g & wrong Q", {
  if (run_test) {
    expect_equal(true_EYgstar - GCOMP_EYgstar, -0.054786)
    expect_equal(true_EYgstar - TMLE_EYgstar,  -0.0005979065)
    expect_equal(tmle_est$estimates[["SE.TMLE"]], 0.0009613481)
  }
})

## -------------------------------------------------
## CORRECT g and CORRECT Q
## 1. sequential GCOMP
## 2. TMLE
## -------------------------------------------------
Qforms <- rep.int("Qkplus1 ~ L1 + L2 + L3 + A + A.tminus1 + L1.tminus1 + L2.tminus1 + L3.tminus1", 2)
# Qmodels <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
lnr_inter <- Lrnr_define_interactions$new(interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
lrn_glm <- Lrnr_glm_fast$new()
Qmodels <- Pipeline$new(lnr_inter, lrn_glm)
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
    expect_equal(true_EYgstar - GCOMP_EYgstar, 0.00130783)
    expect_equal(true_EYgstar - TMLE_EYgstar,  -0.0005873876)
    expect_equal(tmle_est$estimates[["SE.TMLE"]], 0.0009434391)
  }
})

## -------------------------------------------------
## WRONG g and CORRECT Q
## -------------------------------------------------
gform_TRT = "A ~ A.tminus1"
OData <- fitPropensity(OData, gform_TRT = gform_TRT)

IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYgstar <- 1-IPW.St[time == 1, ][, St.directIPW])
cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
## [1] 0.5375674 -- wrong g, IPW is now biased
# IPW bias g^* Astoch:  0.02239061

Qforms <- rep.int("Qkplus1 ~ L1 + L2 + L3 + A + A.tminus1 + L1.tminus1 + L2.tminus1 + L3.tminus1", 2)

# Qmodels <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
lnr_inter <- Lrnr_define_interactions$new(interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
lrn_glm <- Lrnr_glm_fast$new()
Qmodels <- Pipeline$new(lnr_inter, lrn_glm)

tmle_est <- fit_TMLE(OData, tvals = 1, intervened_TRT = "Astoch", Qforms = Qforms, models = Qmodels)
(TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
## [1] 0.5607463 -- TMLE is no longer biased, DR now holds
# TMLE bias g^* Astoch:  -0.0007882782

test_that("Stochastic g^*: IPW is biased and TMLE unbiased for wrong g & correct Q", {
  if (run_test) {
    expect_equal(true_EYgstar - IPW_EYgstar, 0.02239061)
    expect_equal(true_EYgstar - TMLE_EYgstar,  -0.0007882782)
    expect_equal(tmle_est$estimates[["SE.TMLE"]], 0.001030528)
  }
})

# ## -------------------------------------------------
# ## Previous bug (incorrect EIC) was producing wrong results for TMLE under wrong g
# ## WRONG g and CORRECT Q
# ## -------------------------------------------------
# gform_TRT = "A ~ A.tminus1"
# OData <- fitPropensity(OData, gform_TRT = gform_TRT)

# IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
#           directIPW(OData) %$%
#           estimates
# (IPW_EYgstar <- 1-IPW.St[time == 1, ][, St.directIPW])
# cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
# ## [1] 0.5505369 -- wrong g, IPW is now biased
# ## IPW bias g^* Astoch:  0.009421124

# Qforms <- rep.int("Qkplus1 ~ L1 + L2 + L3 + A + A.tminus1 + L1.tminus1 + L2.tminus1 + L3.tminus1", 2)
# params <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))

# tmle_est <- fit_TMLE(OData, tvals = 1, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
# (TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
# cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# ## [1] 0.5495039 -- TMLE is biased, DR no longer holds
# # TMLE bias g^* Astoch:  0.01045413
