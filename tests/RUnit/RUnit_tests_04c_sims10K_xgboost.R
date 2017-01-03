## --------------------------------------------------------------------------------------------------------
## Install data.table (most recent version)
# devtools::install_github('Rdatatable/data.table')
## --------------------------------------------------------------------------------------------------------
## Install stremr
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test xgboost learners
# ---------------------------------------------------------------------------
test.xgboost.10Kdata <- function() {
  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  if (reqxgb) {
    `%+%` <- function(a, b) paste0(a, b)
    options(stremr.verbose = TRUE)
    options(stremr.verbose = FALSE)
    require("data.table")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 1,000 IDs
    Odat_DT <- Odat_DT[ID %in% (1:1000), ]
    setkeyv(Odat_DT, cols = c("ID", "t"))

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

    # ------------------------------------------------------------------
    # Propensity score models for Treatment, Censoring & Monitoring
    # ------------------------------------------------------------------
    gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
    stratify_TRT <- list(
      TI=c("t == 0L",                                            # MODEL TI AT t=0
           "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
           "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
           "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
          ))

    gform_CENS <- c("C ~ highA1c + t")
    gform_MONITOR <- "N ~ 1"

    # ----------------------------------------------------------------
    # IMPORT DATA
    # ----------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost glm
    # ----------------------------------------------------------------
    set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit.method = "cv", fold_column = "fold_ID")

    # set_all_stremr_options(estimator = "xgboost_glm")
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv1 <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv2 <- survNPMSM(wts.St.dhigh, OData)

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    openFile = FALSE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost gbm
    # ----------------------------------------------------------------
    set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm")
    # set_all_stremr_options(estimator = "xgboost_gbm")
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv1 <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv2 <- survNPMSM(wts.St.dhigh, OData)

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    openFile = FALSE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH randomForest (not implemented)
    # ----------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "drf")
    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                         stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

    # wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    # surv1 <- survNPMSM(wts.St.dlow, OData)
    # wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    # surv2 <- survNPMSM(wts.St.dhigh, OData)

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 openFile = FALSE,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit.method = "cv", fold_column = "fold_ID")

    params = list(fit.package = "xgboost", fit.algorithm = "gbm", family = "quasibinomial") # , objective = "reg:logistic"
        # ntrees = 20, learn_rate = 0.1, sample_rate = 0.9, col_sample_rate = 0.9, balance_classes = TRUE)
    # params = list(fit.package = "h2o", fit.algorithm = "randomForest", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)
    t.surv <- c(10)
    Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
    tmle_est$estimates[]

    gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
    gcomp_est$estimates[]

    # ---------------------------------------------------------------------------------------------------------
    # Parallel TMLE w/ xgboost gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    require("doParallel")
    registerDoParallel(cores = 6)
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit.method = "cv", fold_column = "fold_ID")
    params = list(fit.package = "xgboost", fit.algorithm = "gbm", family = "quasibinomial") # , objective = "reg:logistic"
    t.surv <- c(5,6,7,8,9,10)
    Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    tmle_est <- fitTMLE(OData, t_periods = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params,
                        stratifyQ_by_rule = FALSE,
                        parallel = TRUE)
    tmle_est$estimates[]
   est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates   type
    # 1:     TMLE  5 0.1507051 0.8492949           TRUE              0 pooled
    # 2:     TMLE  6 0.1769973 0.8230027           TRUE              0 pooled
    # 3:     TMLE  7 0.2167468 0.7832532           TRUE              0 pooled
    # 4:     TMLE  8 0.2371270 0.7628730           TRUE              0 pooled
    # 5:     TMLE  9 0.2557820 0.7442180           TRUE              0 pooled
    # 6:     TMLE 10 0.2731178 0.7268822           TRUE              0 pooled
    #        TMLE_Var    TMLE_SE rule.name
    # 1: 0.0002323788 0.01524398 gTI.dhigh
    # 2: 0.0002852661 0.01688982 gTI.dhigh
    # 3: 0.0003338560 0.01827173 gTI.dhigh
    # 4: 0.0003632433 0.01905894 gTI.dhigh
    # 5: 0.0003959690 0.01989897 gTI.dhigh
    # 6: 0.0004183633 0.02045393 gTI.dhigh

  }
}

test.xgboost.grid.10Kdata <- function() {
    `%+%` <- function(a, b) paste0(a, b)
    # options(stremr.verbose = TRUE)
    options(stremr.verbose = FALSE)
    require("data.table")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 1,000 IDs
    Odat_DT <- Odat_DT[ID %in% (1:1000), ]
    setkeyv(Odat_DT, cols = c("ID", "t"))

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
    # IMPORT DATA
    # ----------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost gbm
    # ---------------------------------------------------------------------------------------------------------
    set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit.method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(estimator = "xgboost_glm")

    hyper_params = list(eta = c(0.3, 0.1,0.01),
                      max_depth = c(4,6,8,10),
                      max_delta_step = c(0,1),
                      subsample = 1,
                      scale_pos_weight = 1)

    Qgridmodels <- longGriDiSL::defGrid(estimator = "xgboost_gbm",
                           family = "quasibinomial",
                           search_criteria = list(strategy = "RandomDiscrete", max_models = 10),
                           param_grid = hyper_params,
                           early_stopping_rounds = 10,
                           seed = 123456)

    params = list(models = Qgridmodels) # , objective = "reg:logistic"
    # params = list(fit.package = "xgboost", fit.algorithm = "gbm", family = "quasibinomial") # , objective = "reg:logistic"
    # ntrees = 20, learn_rate = 0.1, sample_rate = 0.9, col_sample_rate = 0.9, balance_classes = TRUE)
    # params = list(fit.package = "h2o", fit.algorithm = "randomForest", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)

    t.surv <- c(10)
    Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
    gcomp_est$estimates[]
# >     gcomp_est$estimates[]
#    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates   type
# 1:    GCOMP 10 0.3431796 0.6568204          FALSE             11 pooled
#    rule.name
# 1: gTI.dhigh

    # ---------------------------------------------------------------------------------------------------------
    # Parallel GCOMP w/ xgboost grid search gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    require("doParallel")
    registerDoParallel(cores = 6)
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit.method = "cv", fold_column = "fold_ID")
    t.surv <- c(5,6,7,8,9,10)
    Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE, parallel = TRUE)
    gcomp_est$estimates[]
    #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates   type
    # 1:    GCOMP  5 0.1806779 0.8193221          FALSE              6 pooled
    # 2:    GCOMP  6 0.2171176 0.7828824          FALSE              7 pooled
    # 3:    GCOMP  7 0.2578594 0.7421406          FALSE              8 pooled
    # 4:    GCOMP  8 0.3101151 0.6898849          FALSE              9 pooled
    # 5:    GCOMP  9 0.3300235 0.6699765          FALSE             10 pooled
    # 6:    GCOMP 10 0.3431796 0.6568204          FALSE             11 pooled
    #    rule.name
    # 1: gTI.dhigh
    # 2: gTI.dhigh
    # 3: gTI.dhigh
    # 4: gTI.dhigh
    # 5: gTI.dhigh
    # 6: gTI.dhigh

}



