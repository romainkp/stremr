## --------------------------------------------------------------------------------------------
## Install stremr (LTMLE with stochastic interventions):
# library(devtools)
# install_github("osofr/stremr")
## Install simcausal (simulations with longitudinal data):
# install.packages("simcausal")
## --------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------
## Simulate some data w/ 2 time-points
## Evaluate the truth under stochastic intervention g^* (also under two static interventions)
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
    action("A0", nodes = node("A", distr = "rbern", prob = 0)) +
    action("A1", nodes = node("A", distr = "rbern", prob = 1)) +
    action("Astoch", nodes = node("A", distr = "rbern", prob = 0.3))

  Pstar.Astoch <- data.table(sim(Dset, actions = "Astoch", n = 1000000)[[1]])
  Pstar.A0 <- data.table(sim(Dset, actions = "A0", n = 1000000)[[1]])
  Pstar.A1 <- data.table(sim(Dset, actions = "A1", n = 1000000)[[1]])

  Pstar.Astoch[, lapply(.SD, mean)]
  Pstar.A0[, lapply(.SD, mean)]
  Pstar.A1[, lapply(.SD, mean)]

  (true_EYgstar <- mean(Pstar.Astoch[, Y]))
  ## [1] 0.63725
  (true_EYA0 <- mean(Pstar.A0[, Y]))
  ## [1] 0.591571
  (true_EYA1 <- mean(Pstar.A1[, Y]))
  ## [1] 0.742266

  ## **** UNCOMMENT TO SAVE THIS DATASET ****
  ## WILL SAVE TO CURRENT WORKING DIRECTORY
  # notrun.save.example.data.05 <- function(dt_1t) {
  #   require("tools")
  #   save(dt_1t, compress = TRUE, file = "./dt_1t.rda", compression_level = 9)
  #   resaveRdaFiles("./dt_1t.rda", compress = "bzip2")
  # }

  # attributes(dt_1t)[["true_EYgstar"]] <- true_EYgstar
  # attributes(dt_1t)[["true_EYA0"]] <- true_EYA0
  # attributes(dt_1t)[["true_EYA1"]] <- true_EYA1
  # notrun.save.example.data.05(dt_1t)
}

# library("magrittr")
# library("data.table")

# # ## --------------------------------------------------------------------------------------------
# # ## Perform estimation with LTMLE for stochastic intervention g^*
# # ## Using simulated observed data
# # ## --------------------------------------------------------------------------------------------
# # ## load 2 time-point data example
# # data(dt_2t)
# # true_EYgstar <- attributes(dt_2t)[["true_EYgstar"]]
# # true_EYA0 <- attributes(dt_2t)[["true_EYA0"]]
# # true_EYA1 <- attributes(dt_2t)[["true_EYA1"]]

# # # options(stremr.verbose = TRUE)

# ## have to create a time variable, even though its a constant for all observations
# dt_1t[, "t" := 0]

# ## define the counterfactual A^*
# dt_1t[, ("Astoch") := 0.3][, ("A0") := 0][, ("A1") := 1]

# ## add some lags
# OData <- importData(dt_1t, ID = "ID", t = "t", covars = c("L1", "L2", "L3", "L4", "L5", "L6"), TRT = "A", OUTCOME = "Y")

# ## -------------------------------------------------
# ## IPW with correct g
# ## -------------------------------------------------
# gform_TRT = "A ~ L1 + L2 + L3 + L4 + L5 + L6"
# ## If you want to stratify the propensity score fit by some baseline covariate, use this option.
# ## Otherwise just leave it out and the propensity score will use all subjects to fit a
# ## single propensity score model
# stratify_TRT <- list(A=c("L1 == 0", "L1 == 1"))
# ## Define the library of candidate estimators for the propensity score model
# models_TRT <- defModel(estimator = "speedglm__glm")
# OData <- fitPropensity(OData, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, models_TRT = models_TRT)

# IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
#           directIPW(OData) %$%
#           estimates
# (IPW_EYgstar <- 1-IPW.St[time == 0, ][, St.directIPW])
# cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
# ## [1]
# # IPW bias g^* Astoch:  -3.935335e-05

# IPW.St <- getIPWeights(OData, intervened_TRT = "A0") %>%
#           directIPW(OData) %$%
#           estimates
# (IPW_EYA0 <- 1-IPW.St[time == 0, ][, St.directIPW])
# cat("\nIPW bias A0: ", true_EYA0-IPW_EYA0, "\n")
# ## [1]
# ## IPW bias A0:

# IPW.St <- getIPWeights(OData, intervened_TRT = "A1") %>%
#           directIPW(OData) %$%
#           estimates
# (IPW_EYA1 <- 1-IPW.St[time == 0, ][, St.directIPW])
# cat("\nIPW bias A1: ", true_EYA1-IPW_EYA1, "\n")
# ## [1]
# ## IPW bias A1:

# ## -------------------------------------------------
# ## CORRECT g and CORRECT Q
# ## 1. sequential GCOMP; 2. TMLE
# ## -------------------------------------------------
# Qforms <- "Qkplus1 ~ A + L1 + L2 + L3 + L4 + L5 + L6"
# ## UNCOMMENT TO run with xgboost GBM
# # params <- gridisl::defModel(estimator = "xgboost__gbm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
# params <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
# gcomp_est <- fit_GCOMP(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
# (GCOMP_EYgstar <- 1 - gcomp_est$estimates[, St.GCOMP])
# cat("\nGCOMP bias g^* Astoch: ", true_EYgstar-GCOMP_EYgstar, "\n")
# ## (glm) [1]  -- unbiased
# ## GCOMP bias g^* Astoch:

# tmle_est <- fit_TMLE(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
# (TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
# cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# ## [1]  -- unbiased
# ## TMLE bias g^* Astoch:

# ## -------------------------------------------------
# ## CORRECT g and WRONG Q
# ## 1. GCOMP for stochastic intervention; 2. TMLE
# ## -------------------------------------------------
# Qforms <- "Qkplus1 ~ L1 + L2"
# params <- gridisl::defModel(estimator = "speedglm__glm")
# gcomp_est <- fit_GCOMP(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
# (GCOMP_EYgstar <- 1 - gcomp_est$estimates[, St.GCOMP])
# cat("\nGCOMP bias g^* Astoch: ", true_EYgstar-GCOMP_EYgstar, "\n")
# ## [1]  -- very biased
# ## GCOMP bias g^* Astoch:

# tmle_est <- fit_TMLE(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
# (TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
# cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# ##  -- unbiased, same as IPW
# ## TMLE bias g^* Astoch:

# ## -------------------------------------------------
# ## WRONG g and CORRECT Q
# ## -------------------------------------------------
# gform_TRT = "A ~ 1"
# OData <- fitPropensity(OData, gform_TRT = gform_TRT)

# IPW.St <- getIPWeights(OData, intervened_TRT = "Astoch") %>%
#           directIPW(OData) %$%
#           estimates
# (IPW_EYgstar <- 1-IPW.St[time == 0, ][, St.directIPW])
# cat("\nIPW bias g^* Astoch: ", true_EYgstar-IPW_EYgstar, "\n")
# ## [1]  -- wrong g, IPW is now biased
# # IPW bias g^* Astoch:

# Qforms <- "Qkplus1 ~ L1 + L2 + L3 + A "
# params <- gridisl::defModel(estimator = "speedglm__glm", interactions = list(c("A", "L1"), c("A", "L2"), c("A", "L2")))
# tmle_est <- fit_TMLE(OData, tvals = 0, intervened_TRT = "Astoch", Qforms = Qforms, models = params)
# (TMLE_EYgstar <- 1 - tmle_est$estimates[, St.TMLE])
# cat("\nTMLE bias g^* Astoch: ", true_EYgstar-TMLE_EYgstar, "\n")
# ## [1]  -- TMLE is no longer biased, DR now holds
# # TMLE bias g^* Astoch:

