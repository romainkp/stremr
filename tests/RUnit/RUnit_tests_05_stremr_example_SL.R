# devtools::install_github("Rdatatable/data.table")
# --------------------------------------------------------------------------------------------------------
# Install h2o (most recent version)
# --------------------------------------------------------------------------------------------------------
# if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
# if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
# # Next, we download packages that H2O depends on.
# if (! ("methods" %in% rownames(installed.packages()))) { install.packages("methods") }
# if (! ("statmod" %in% rownames(installed.packages()))) { install.packages("statmod") }
# if (! ("stats" %in% rownames(installed.packages()))) { install.packages("stats") }
# if (! ("graphics" %in% rownames(installed.packages()))) { install.packages("graphics") }
# if (! ("RCurl" %in% rownames(installed.packages()))) { install.packages("RCurl") }
# if (! ("jsonlite" %in% rownames(installed.packages()))) { install.packages("jsonlite") }
# if (! ("tools" %in% rownames(installed.packages()))) { install.packages("tools") }
# if (! ("utils" %in% rownames(installed.packages()))) { install.packages("utils") }
# # Now we download, install and initialize the H2O package for R.
# install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/rel-turchin/9/R")))
# --------------------------------------------------------------------------------------------------------
# Install h2oEnsemble (most recent stable version 1.8)
# --------------------------------------------------------------------------------------------------------
# install.packages("https://h2o-release.s3.amazonaws.com/h2o-ensemble/R/h2oEnsemble_0.1.8.tar.gz", repos = NULL)
#   Install h2oEnsemble (dev version):
# devtools::install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")

test.h2oEnsemble <- function() {
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
    # This is set to run fairly quickly, increase max_runtime_secs
    # or max_models to cover more of the hyperparameter space.
    # Also, you can expand the hyperparameter space of each of the
    # algorithms by modifying the hyper param code below.

    # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 20)
    # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 120)
    # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 5*120)
    # search_criteria <- list(strategy = "RandomDiscrete", max_models = 42, max_runtime_secs = 28800)
    # search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "AUTO", stopping_tolerance = 0.001, stopping_rounds = 10)
    # search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "misclassification", stopping_tolerance = 0.00001, stopping_rounds = 5)
    # nfolds <- 5
    # nfolds <- 10

    ntrees_opt <- c(100, 200, 300, 500)
    mtries_opt <- 8:20
    max_depth_opt <- c(5, 10, 15, 20, 25)
    sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
    col_sample_rate_per_tree_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    balance_classes_opt <- c(TRUE, FALSE)
    RF_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 20),
                            ntrees = ntrees_opt,
                            # mtries = mtries_opt,
                            max_depth = max_depth_opt,
                            sample_rate = sample_rate_opt,
                            col_sample_rate_per_tree = col_sample_rate_per_tree_opt,
                            balance_classes = balance_classes_opt)

    ntrees_opt <- c(100, 200, 300, 500)
    learn_rate_opt <- c(0.005, 0.01, 0.03, 0.06)
    max_depth_opt <- c(3, 4, 5, 6, 9)
    sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
    col_sample_rate_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    balance_classes_opt <- c(TRUE, FALSE)
    GBM_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 20),
                             ntrees = ntrees_opt,
                             learn_rate = learn_rate_opt,
                             max_depth = max_depth_opt,
                             sample_rate = sample_rate_opt,
                             col_sample_rate = col_sample_rate_opt,
                             balance_classes = balance_classes_opt)

    # activation_opt <- c("Rectifier", "RectifierWithDropout", "Maxout", "MaxoutWithDropout")
    # hidden_opt <- list(c(10,10), c(20,15), c(50,50,50))
    # l1_opt <- c(0, 1e-3, 1e-5)
    # l2_opt <- c(0, 1e-3, 1e-5)
    # DL_hyper_params <- list(activation = activation_opt,
    #                               hidden = hidden_opt,
    #                               l1 = l1_opt,
    #                               l2 = l2_opt)

    h2o.randomForest.1 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
    h2o.randomForest.2 <- function(..., ntrees = 200, sample_rate = 0.75, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
    h2o.randomForest.3 <- function(..., ntrees = 200, sample_rate = 0.85, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
    h2o.gbm.1 <- function(..., ntrees = 100, nbins = 100, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
    h2o.gbm.2 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
    h2o.gbm.3 <- function(..., ntrees = 100, max_depth = 10, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
    h2o.gbm.4 <- function(..., ntrees = 100, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
    h2o.gbm.5 <- function(..., ntrees = 200, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
    h2o.gbm.6 <- function(..., ntrees = 200, col_sample_rate = 0.7, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
    # h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
    # h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
    # h2o.deeplearning.3 <- function(..., hidden = c(500,500), activation = "RectifierWithDropout", seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
    # h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", seed = 1) h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
    # h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", seed = 1) h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
    # learner <- c("h2o.randomForest.1", "h2o.deeplearning.1", "h2o.deeplearning.2")

    h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
    h2o.glm.2 <- function(..., x = "highA1c", alpha = 0.0) h2o.glm.wrapper(..., x = x, alpha = alpha)
    h2o.glm.3 <- function(..., alpha = 1.0) h2o.glm.wrapper(..., alpha = alpha)
    # learner <- c("h2o.glm.1", "h2o.glm.2", "h2o.glm.3")
                 # "h2o.randomForest.1", "h2o.randomForest.2", "h2o.randomForest.3",
                 # "h2o.gbm.1", "h2o.gbm.2", "h2o.gbm.3", "h2o.gbm.4", "h2o.gbm.5", "h2o.gbm.6")
                # "h2o.deeplearning.1", "h2o.deeplearning.2", "h2o.deeplearning.3"
    # metalearner <- "h2o.glm_nn"
    # family <- "binomial"
    # h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)

    glm_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
                             alpha = c(0,1,seq(0.1,0.9,0.1)), lambda = c(0,1e-7,1e-5,1e-3,1e-1))

    SLparams = list( # search_criteria = list(strategy = "RandomDiscrete", max_runtime_secs = 20),
                     grid.algorithm = c("glm"),
                     # grid.algorithm = c("glm", "randomForest"),
                     # learner = c("h2o.glm.2"),
                     learner = c("h2o.glm.wrapper"),
                     # algorithm = c("glm", "randomForest", "gbm", "deeplearning"),
                     metalearner = "h2o.glm_nn",
                     nfolds = 2,
                     # nfolds = 5,
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
    make_report_rmd(OData, openFile = FALSE,
                    # AddFUPtables = TRUE,
                    # RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = FALSE),
                    # WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)

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
    tmle_est[]

    h2o::h2o.shutdown(prompt = FALSE)
}