## --------------------------------------------------------------------------------------------
## Install stremr (LTMLE with stochastic interventions):
# library(devtools)
# install_github("osofr/stremr")
## Install simcausal (simulations with longitudinal data):
# install.packages("simcausal")
## --------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------
## Simulate some data w/ 1 time-point
## Evaluate the truth under two static interventions and one stochastic intervention g^*
## --------------------------------------------------------------------------------------------
sim_stochastic_test_data_1t <- function() {
  library("data.table")
  library("simcausal")
  set.seed(123456)

  D <- DAG.empty()
  D <- D +
    node("L1", distr = "rbern", prob = 0.5) +
    node("L2", distr = "rbern", prob = 0.5) +
    node("L3", distr = "rbern", prob = 0.5) +
    node("L4", distr = "rbern", prob = plogis(-0.5 + 0.5*L1 + 0.5*L2)) +
    node("L5", distr = "rbern", prob = plogis(-0.5 + 0.5*L1 + 0.5*L2)) +
    node("L6", distr = "rbern", prob = plogis(-0.5 + 0.5*L1 + 0.5*L3)) +
    node("A",  distr = "rbern", prob = plogis(0.5*L1 + 0.5*L2 + 0.5*L3)) +
    node("Y",  distr = "rbern", prob = plogis(0.5*A + 0.5*L2 + 0.5*L3 + 0.25*L4 + 0.25*L5 + 0.25*L6))

  Dset <- set.DAG(D)
  dt_1t <- data.table(sim(Dset, n = 50000))
  dt_1t[, lapply(.SD, mean)]

  Dset <- Dset +
    ## static intervention that sets everyone's treatment to 0:
    action("A0", nodes = node("A", distr = "rbern", prob = 0)) +
    ## static intervention that sets everyone's treatment to 1:
    action("A1", nodes = node("A", distr = "rbern", prob = 1)) +
    ## stochastic intervention where P(A^*=1)=0.3:
    action("Astoch", nodes = node("A", distr = "rbern", prob = 0.3))

  Pstar.Astoch <- data.table(sim(Dset, actions = "Astoch", n = 1000000)[[1]])
  Pstar.A0 <- data.table(sim(Dset, actions = "A0", n = 1000000)[[1]])
  Pstar.A1 <- data.table(sim(Dset, actions = "A1", n = 1000000)[[1]])

  Pstar.Astoch[, lapply(.SD, mean)]
  Pstar.A0[, lapply(.SD, mean)]
  Pstar.A1[, lapply(.SD, mean)]

  print(true_EYgstar <- mean(Pstar.Astoch[, Y]))
  ## [1] 0.724942
  print(true_EYA0 <- mean(Pstar.A0[, Y]))
  ## [1] 0.697425
  print(true_EYA1 <- mean(Pstar.A1[, Y]))
  ## [1] 0.789542

  ## **** UNCOMMENT TO SAVE THIS DATASET ****
  ## WILL SAVE TO CURRENT WORKING DIRECTORY
  notrun.save.example.data.05 <- function(dt_1t) {
    require("tools")
    save(dt_1t, compress = TRUE, file = "./dt_1t.rda", compression_level = 9)
    resaveRdaFiles("./dt_1t.rda", compress = "bzip2")
  }

  attributes(dt_1t)[["true_EYgstar"]] <- true_EYgstar
  attributes(dt_1t)[["true_EYA0"]] <- true_EYA0
  attributes(dt_1t)[["true_EYA1"]] <- true_EYA1
  # notrun.save.example.data.05(dt_1t)
  return(dt_1t)
}

library("stremr")
library("magrittr")
library("data.table")

## --------------------------------------------------------------------------------------------
## obtain the observed data with single time-point (and truth(s) saved as attributes)
## --------------------------------------------------------------------------------------------
# dt_1t <- sim_stochastic_test_data_1t()
# dt_1t <- copy(dt_1t)
data(dt_1t)
## --------------------------------------------------------------------------------------------
## Perform estimation with stremr (LTMLE) for stochastic intervention g^*
## Using simulated observed data
## --------------------------------------------------------------------------------------------
true_EYgstar <- attributes(dt_1t)[["true_EYgstar"]]
true_EYA0 <- attributes(dt_1t)[["true_EYA0"]]
true_EYA1 <- attributes(dt_1t)[["true_EYA1"]]

options(stremr.verbose = FALSE)
options(gridisl.verbose = FALSE)
options(sl3.verbose = FALSE)
options(condensier.verbose = FALSE)
# options(stremr.verbose = TRUE)
# options(gridisl.verbose = TRUE)
# options(sl3.verbose = TRUE)
# options(condensier.verbose = TRUE)

## have to create a time variable, even though its a constant for all observations
dt_1t[, "t" := 0]

## define the counterfactual A^* (stochastic, static A=0, static A=1)
dt_1t[, ("Astoch") := 0.3][, ("A0") := 0][, ("A1") := 1]

## add some lags
OData <- importData(dt_1t, ID = "ID", t = "t", covars = c("L1", "L2", "L3", "L4", "L5", "L6"), TRT = "A", OUTCOME = "Y")

## -------------------------------------------------
## specify correct g and fit TRT model / propensity score
## -------------------------------------------------
gform_TRT = "A ~ L1 + L2 + L3 + L4 + L5 + L6"

## If you want to stratify the propensity score fit by some baseline covariate, use this option.
## Otherwise just leave it out and the propensity score will use all subjects to fit a
## single propensity score model
stratify_TRT <- list(A=c("L1 == 0", "L1 == 1"))
## Define the library of candidate estimators for the propensity score model
models_TRT <- defModel(estimator = "speedglm__glm")
OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = models_TRT)

## -------------------------------------------------
## IPW with correct g for stochastic intervention "Astoch"
## -------------------------------------------------
IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYgstar <- 1-IPW.St[time == 0, ][, St.directIPW])
cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
# IPW bias g^* Astoch: 0.000520762

## -------------------------------------------------
## IPW with correct g for static intervention "A0"
## -------------------------------------------------
IPW.St <- getIPWeights(OData, intervened_TRT = "A0") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYA0 <- 1-IPW.St[time == 0, ][, St.directIPW])
cat("\nIPW bias A0: ", true_EYA0-IPW_EYA0, "\n")
# IPW bias A0:  -0.001090636

## -------------------------------------------------
## IPW with correct g for static intervention "A1"
## -------------------------------------------------
IPW.St <- getIPWeights(OData, intervened_TRT = "A1") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYA1 <- 1-IPW.St[time == 0, ][, St.directIPW])
cat("\nIPW bias A1: ", true_EYA1-IPW_EYA1, "\n")
# IPW bias A1:  0.004695718 

## -------------------------------------------------
## CORRECT g and CORRECT Q for estimators for stochastic intervention:
## 1. sequential GCOMP; 2. TMLE
## -------------------------------------------------
## specify the outcome model
Qforms <- "Qkplus1 ~ A + L1 + L2 + L3 + L4 + L5 + L6"
## UNCOMMENT TO run with xgboost GBM
# params <- gridisl::defModel(estimator = "xgboost__gbm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
params <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
gcomp_est <- fit_GCOMP(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
(GCOMP_EYgstar <- 1 - gcomp_est$estimates[, St.GCOMP])
cat("\nGCOMP bias g^* Astoch: ", true_EYgstar-GCOMP_EYgstar, "\n")
# GCOMP bias g^* Astoch:  -0.0001249876 

tmle_est <- fit_TMLE(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
(TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# TMLE bias g^* Astoch:  0.0003941923 

## -------------------------------------------------
## CORRECT g and WRONG Q for estimators for stochastic interventions:
## 1. sequential GCOMP; 2. TMLE
## -------------------------------------------------
Qforms <- "Qkplus1 ~ L1 + L2"
params <- gridisl::defModel(estimator = "speedglm__glm")
gcomp_est <- fit_GCOMP(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
(GCOMP_EYgstar <- 1 - gcomp_est$estimates[, St.GCOMP])
cat("\nGCOMP bias g^* Astoch: ", true_EYgstar-GCOMP_EYgstar, "\n")
## [1]  -- very biased
# GCOMP bias g^* Astoch:  -0.030698 

tmle_est <- fit_TMLE(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
(TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# TMLE bias g^* Astoch:  0.0004701075 

## -------------------------------------------------
## WRONG g and CORRECT Q
## -------------------------------------------------
gform_TRT = "A ~ 1"
OData <- fitPropensity(OData, gform_TRT = gform_TRT)

IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
          directIPW(OData) %$%
          estimates
(IPW_EYgstar <- 1-IPW.St[time == 0, ][, St.directIPW])
cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
## [1]  -- wrong g, IPW is now biased
# IPW bias g^* Astoch:  0.01100976 

Qforms <- "Qkplus1 ~ L1 + L2 + L3 + A "
params <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
tmle_est <- fit_TMLE(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
(TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# TMLE bias g^* Astoch:  -8.566571e-05 

