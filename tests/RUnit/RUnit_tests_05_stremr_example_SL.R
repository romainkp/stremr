# --------------------------------------------------------------------------------------------------------
# Install data.table (most recent version)
# --------------------------------------------------------------------------------------------------------
# devtools::install_github('Rdatatable/data.table')
# --------------------------------------------------------------------------------------------------------
# Install h2o (most recent version)
# --------------------------------------------------------------------------------------------------------
# if ("package:h2o" %in% search()) detach("package:h2o", unload=TRUE)
# if ("h2o" %in% rownames(installed.packages())) remove.packages("h2o")
# # Next, download H2O package dependencies:
# pkgs <- c("methods","statmod","stats","graphics","RCurl","jsonlite","tools","utils")
# new.pkgs <- setdiff(pkgs, rownames(installed.packages()))
# if (length(new.pkgs)) install.packages(new.pkgs)
# # Download and install the H2O package for R:
# install.packages("h2o", type="source", repos=(c("https://s3.amazonaws.com/h2o-release/h2o/master/3636/R")))
# --------------------------------------------------------------------------------------------------------
# Install h2oEnsemble (most recent stable version 1.8)
# --------------------------------------------------------------------------------------------------------
# install.packages("https://h2o-release.s3.amazonaws.com/h2o-ensemble/R/h2oEnsemble_0.1.8.tar.gz", repos = NULL)
#   Install h2oEnsemble (dev version):
# devtools::install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")
# --------------------------------------------------------------------------------------------------------
# Install stremr
# --------------------------------------------------------------------------------------------------------
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)

test.h2oEnsemble <- function() {
    reqh2o <- requireNamespace("h2o", quietly = TRUE)
    reqSL <- requireNamespace("h2oEnsemble", quietly = TRUE)
    if (reqh2o && reqSL) {
        options(stremr.verbose = TRUE)
        `%+%` <- function(a, b) paste0(a, b)
        require("h2o")
        require('h2oEnsemble')
        require("data.table")
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
        # IMPORT DATA / INIT h2o
        # ----------------------------------------------------------------
        set_all_stremr_options(fit.package = "h2o", fit.algorithm = "SuperLearner")
        # set_all_stremr_options(fit.package = "h2o", fit.algorithm = "RF")
        # h2o::h2o.init(nthreads = -1)
        h2o::h2o.init(nthreads = 1)
        # h2o::h2o.init(nthreads = -1, startH2O = FALSE)
        # h2o::h2o.shutdown(prompt = FALSE)
        OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
        # OData$define_CVfolds(nfolds = 10, seed = 23)
        # OData$dat.sVar

        # ------------------------------------------------------------------
        # Fit propensity scores for Treatment, Censoring & Monitoring
        # ------------------------------------------------------------------
        # gform_TRT <- c("TI ~ CVD + highA1c", "TI ~ CVD + highA1c", "TI ~ CVD + highA1c", "TI ~ 1")
        gform_TRT <- c("TI ~ CVD + highA1c")
        stratify_TRT <- list(
          TI=c("t == 0L",                                            # MODEL TI AT t=0
               "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
               "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
               "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
              ))

        gform_CENS <- c("C ~ highA1c + t")
        gform_MONITOR <- "N ~ 1"

        # Random Grid Search (e.g. 120 second maximum)
        # This is set to run fairly quickly, increase max_runtime_secs or max_models to cover more of the hyperparameter space.
        # You can expand the hyperparameter space of each of the algorithms by modifying the hyper param code below.
        # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 20)
        # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 120)
        # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 5*120)
        # search_criteria <- list(strategy = "RandomDiscrete", max_models = 42, max_runtime_secs = 28800)
        # search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "AUTO", stopping_tolerance = 0.001, stopping_rounds = 10)
        # search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "misclassification", stopping_tolerance = 0.00001, stopping_rounds = 5)

        glm_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
                                 alpha = c(0,1,seq(0.1,0.9,0.1)),
                                 lambda = c(0,1e-7,1e-5,1e-3,1e-1))

        RF_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 20),
                                ntrees = c(100, 200, 300, 500),
                                mtries = 8:20,
                                max_depth = c(5, 10, 15, 20, 25),
                                sample_rate = c(0.7, 0.8, 0.9, 1.0),
                                col_sample_rate_per_tree = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                                balance_classes = c(TRUE, FALSE))

        GBM_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 20),
                                 ntrees = c(100, 200, 300, 500),
                                 learn_rate = c(0.005, 0.01, 0.03, 0.06),
                                 max_depth = c(3, 4, 5, 6, 9),
                                 sample_rate = c(0.7, 0.8, 0.9, 1.0),
                                 col_sample_rate = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                                 balance_classes = c(TRUE, FALSE))

        DL_hyper_params <- list(activation = c("Rectifier", "RectifierWithDropout", "Maxout", "MaxoutWithDropout"),
                                hidden = list(c(10,10), c(20,15), c(50,50,50)),
                                l1 = c(0, 1e-3, 1e-5),
                                l2 = c(0, 1e-3, 1e-5))

        h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
        h2o.glm.2 <- function(..., x = "highA1c", alpha = 0.0) h2o.glm.wrapper(..., x = x, alpha = alpha)
        h2o.glm.3 <- function(..., alpha = 1.0) h2o.glm.wrapper(..., alpha = alpha)

        h2o.randomForest.1 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
        h2o.randomForest.2 <- function(..., ntrees = 200, sample_rate = 0.75, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
        h2o.randomForest.3 <- function(..., ntrees = 200, sample_rate = 0.85, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)

        h2o.gbm.1 <- function(..., ntrees = 100, nbins = 100, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
        h2o.gbm.2 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
        h2o.gbm.3 <- function(..., ntrees = 100, max_depth = 10, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
        h2o.gbm.4 <- function(..., ntrees = 100, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
        h2o.gbm.5 <- function(..., ntrees = 200, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
        h2o.gbm.6 <- function(..., ntrees = 200, col_sample_rate = 0.7, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)

        h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
        h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
        h2o.deeplearning.3 <- function(..., hidden = c(500,500), activation = "RectifierWithDropout", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)

        # learner <- c("h2o.glm.1", "h2o.randomForest.1", "h2o.gbm.1", "h2o.deeplearning.1")
        # learner <- c("h2o.glm.1", "h2o.glm.2", "h2o.glm.3",
        #              "h2o.randomForest.1", "h2o.randomForest.2", "h2o.randomForest.3",
        #              "h2o.gbm.1", "h2o.gbm.2", "h2o.gbm.3", "h2o.gbm.4", "h2o.gbm.5", "h2o.gbm.6",
        #              "h2o.deeplearning.1", "h2o.deeplearning.2", "h2o.deeplearning.3")

        # metalearner <- "h2o.glm_nn"
        # family <- "binomial"
        # h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)

        SLparams = list( # search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 20),
                         grid.algorithm = c("glm"),
                         # grid.algorithm = c("glm", "randomForest"),
                         # learner = c("h2o.glm.2"),
                         learner = c("h2o.glm.wrapper"),
                         metalearner = "h2o.glm_nn",
                         nfolds = 2,
                         # nfolds = 5,
                         family = "binomial",
                         seed = 23,
                         glm = glm_hyper_params
                         # randomForest = RF_hyper_params
                         # gbm = GBM_hyper_params,
                         # deeplearning = DL_hyper_params
                         )

        params_CENS = c(SLparams, save.ensemble = FALSE, ensemble.dir.path = "./h2o-ensemble-model-CENS")
        params_TRT = c(SLparams, save.ensemble = FALSE, ensemble.dir.path = "./h2o-ensemble-model-TRT")
        params_MONITOR = list(fit.package = "speedglm", fit.algorithm = "glm")

        OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                                stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                                params_CENS = params_CENS, params_TRT = params_TRT, params_MONITOR = params_MONITOR)

        # stop("SL end")

        require("magrittr")
        St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly") %>%
                   survNPMSM(OData)  %$%
                   IPW_estimates
        St.dlow

        St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
                    survNPMSM(OData) %$%
                    IPW_estimates
        St.dhigh

        # report.path <- "/set/your/report/path/"
        # file.path = report.path
        if (rmarkdown::pandoc_available(version = "1.12.3"))
            make_report_rmd(OData, openFile = FALSE,
                        # AddFUPtables = TRUE,
                        # MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = FALSE),
                        # WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                        file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", y_legend = 0.95)

        # ---------------------------------------------------------------------------------------------------------
        # ERROR CHECK: no hyper params for randomForest:
        # ---------------------------------------------------------------------------------------------------------
        # options(stremr.verbose = TRUE)
        set_all_stremr_options(fit.package = "h2o", fit.algorithm = "SuperLearner")
        SLparams_testERROR = list(
                         grid.algorithm = c("glm", "randomForest"),
                         nfolds = 3,
                         glm = glm_hyper_params)
        params_CENS = SLparams_testERROR
        params_TRT = list(fit.package = "speedglm", fit.algorithm = "glm")
        params_MONITOR = list(fit.package = "speedglm", fit.algorithm = "glm")
        OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                                stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                                params_CENS = params_CENS, params_TRT = params_TRT, params_MONITOR = params_MONITOR)

        # ---------------------------------------------------------------------------------------------------------
        # TMLE / GCOMP
        # ---------------------------------------------------------------------------------------------------------
        t.surv <- c(5)
        Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
        params = list(fit.package = "speedglm", fit.algorithm = "glm")
        # params = list(fit.package = "h2o", fit.algorithm = "RF", ntrees = 100,
        #               learn_rate = 0.05, sample_rate = 0.8,
        #               col_sample_rate = 0.8, balance_classes = TRUE)
        Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
        tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
        tmle_est$estimates[]

        h2o::h2o.shutdown(prompt = FALSE)
    }
}