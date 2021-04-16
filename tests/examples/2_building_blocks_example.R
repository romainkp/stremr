options(stremr.verbose = TRUE)
require("data.table")

# ----------------------------------------------------------------------
# Simulated Data
# ----------------------------------------------------------------------
data(OdataNoCENS)
OdataDT <- as.data.table(OdataNoCENS, key=c("ID", "t"))

# define lagged N, first value is always 1 (always monitored at the first time point):
OdataDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
OdataDT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]

# ----------------------------------------------------------------------
# Define intervention (always treated):
# ----------------------------------------------------------------------
OdataDT[, ("TI.set1") := 1L]
OdataDT[, ("TI.set0") := 0L]

# ----------------------------------------------------------------------
# Import Data
# ----------------------------------------------------------------------
OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "N.tminus1"),
                    CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")

# ----------------------------------------------------------------------
# Look at the input data object
# ----------------------------------------------------------------------
print(OData)

# ----------------------------------------------------------------------
# Access the input data
# ----------------------------------------------------------------------
get_data(OData)

# ----------------------------------------------------------------------
# Model the Propensity Scores
# ----------------------------------------------------------------------
gform_CENS <- "C ~ highA1c + lastNat1"
gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
gform_MONITOR <- "N ~ 1"
stratify_CENS <- list(C=c("t < 16", "t == 16"))

# ----------------------------------------------------------------------
# Fit Propensity Scores
# ----------------------------------------------------------------------
OData <- fitPropensity(OData, gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        gform_MONITOR = gform_MONITOR,
                        stratify_CENS = stratify_CENS)

# ----------------------------------------------------------------------
# IPW Ajusted KM or Saturated MSM
# ----------------------------------------------------------------------
require("magrittr")
AKME.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
             survNPMSM(OData) %$%
             estimates
AKME.St.1

# ----------------------------------------------------------------------
# Bounded IPW
# ----------------------------------------------------------------------
IPW.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
            directIPW(OData)
IPW.St.1[]

# ----------------------------------------------------------------------
# IPW-MSM for hazard
# ----------------------------------------------------------------------
wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, tbreaks = c(1:8,12,16)-1,)
survMSM_res$St

# ----------------------------------------------------------------------
# Sequential G-COMP
# ----------------------------------------------------------------------
t.surv <- c(0:10)
Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
params <- gridisl::defModel(estimator = "speedglm__glm")

\dontrun{
gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "TI.set1",
                          Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)
gcomp_est[]
}
# ----------------------------------------------------------------------
# TMLE
# ----------------------------------------------------------------------
\dontrun{
tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "TI.set1",
                    Qforms = Qforms, models = params, stratifyQ_by_rule = TRUE)
tmle_est[]
}

# ----------------------------------------------------------------------
# Running IPW-Adjusted KM with optional user-specified weights:
# ----------------------------------------------------------------------
addedWts_DT <- OdataDT[, c("ID", "t"), with = FALSE]
addedWts_DT[, new.wts := sample.int(10, nrow(OdataDT), replace = TRUE)/10]
survNP_res_addedWts <- survNPMSM(wts.DT.1, OData, weights = addedWts_DT)

# ----------------------------------------------------------------------
# Multivariate Propensity Score Regressions
# ----------------------------------------------------------------------
gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                        gform_MONITOR = gform_MONITOR)

# ----------------------------------------------------------------------
# Fitting treatment model with Gradient Boosting machines:
# ----------------------------------------------------------------------
\dontrun{
require("h2o")
h2o::h2o.init(nthreads = -1)
gform_CENS <- "C ~ highA1c + lastNat1"
models_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "gbm")
OData <- fitPropensity(OData, gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        models_TRT = models_TRT,
                        gform_MONITOR = gform_MONITOR,
                        stratify_CENS = stratify_CENS)

# Use `H2O-3` distributed implementation of GLM for treatment model estimator:
models_TRT <- sl3::Lrnr_h2o_glm$new(family = "binomial")
OData <- fitPropensity(OData, gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        models_TRT = models_TRT,
                        gform_MONITOR = gform_MONITOR,
                        stratify_CENS = stratify_CENS)

# Use Deep Neural Nets:
models_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "deeplearning")
OData <- fitPropensity(OData, gform_CENS = gform_CENS,
                        gform_TRT = gform_TRT,
                        models_TRT = models_TRT,
                        gform_MONITOR = gform_MONITOR,
                        stratify_CENS = stratify_CENS)
}

# ----------------------------------------------------------------------
# Fitting different models with different algorithms
# Fine tuning modeling with optional tuning parameters.
# ----------------------------------------------------------------------
\dontrun{
params_TRT <- sl3::Lrnr_h2o_grid$new(algorithm = "gbm",
                              ntrees = 50,
                              learn_rate = 0.05,
                              sample_rate = 0.8,
                              col_sample_rate = 0.8,
                              balance_classes = TRUE)
params_CENS <- sl3::Lrnr_glm_fast$new()
params_MONITOR <- sl3::Lrnr_glm_fast$new()
OData <- fitPropensity(OData,
            gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, params_CENS = params_CENS,
            gform_TRT = gform_TRT, params_TRT = params_TRT,
            gform_MONITOR = gform_MONITOR, params_MONITOR = params_MONITOR)
}

# ----------------------------------------------------------------------
# Running TMLE based on the previous fit of the propensity scores.
# Also applying Random Forest to estimate the sequential outcome model
# ----------------------------------------------------------------------
\dontrun{
t.surv <- c(0:5)
Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
models <- sl3::Lrnr_h2o_grid$new(algorithm = "randomForest",
                           ntrees = 100, learn_rate = 0.05, sample_rate = 0.8,
                           col_sample_rate = 0.8, balance_classes = TRUE)
tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "TI.set1",
            Qforms = Qforms, models = models,
            stratifyQ_by_rule = TRUE)
}

\dontrun{
t.surv <- c(0:5)
Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
models <- sl3::Lrnr_h2o_grid$new(algorithm = "randomForest",
                           ntrees = 100, learn_rate = 0.05, sample_rate = 0.8,
                           col_sample_rate = 0.8, balance_classes = TRUE)
tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "TI.set1",
            Qforms = Qforms, models = models,
            stratifyQ_by_rule = FALSE)
}

