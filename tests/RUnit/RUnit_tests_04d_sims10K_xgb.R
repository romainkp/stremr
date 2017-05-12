## --------------------------------------------------------------------------------------------------------
## Install data.table (most recent version)
# devtools::install_github('Rdatatable/data.table')
## --------------------------------------------------------------------------------------------------------
## Install stremr
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test reporting with a mix of h2o / xgb models, using CV
# ---------------------------------------------------------------------------
test.GRID.h2o.xgboost.10Kdata <- function() {
  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  reqh2o <- requireNamespace("h2o", quietly = TRUE)
  if (reqxgb && reqh2o) {
    `%+%` <- function(a, b) paste0(a, b)
    # library("stremr")
    # options(stremr.verbose = TRUE)
    options(stremr.verbose = FALSE)
    options(gridisl.verbose = TRUE)
    # options(gridisl.verbose = FALSE)
    library("data.table")
    library("h2o")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 1,000 IDs
    # Odat_DT <- Odat_DT[ID %in% (1:1000), ]
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
    OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    params_g <- gridisl::defModel(estimator = "xgboost__glm",
                                  family = "binomial",
                                  nrounds = 10,
                                  early_stopping_rounds = 2)

    OData <- fitPropensity(OData,
                            gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            models_CENS = params_g, models_TRT = params_g, models_MONITOR = params_g,
                            fit_method = "cv",
                            fold_column = "fold_ID")

    ## ----------------------------------------------------------------------------------------------------
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.periods.RDs)+1))
    # set_all_stremr_options(estimator = "speedglm__glm")
    tmle_est_dlow <- fitTMLE(OData, tvals = t.periods.RDs, intervened_TRT = "gTI.dlow",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        estimator = "speedglm__glm",
                        family = "quasibinomial")

    tmle_est_dlow[["estimates"]]

    tmle_est_dhigh <- fitTMLE(OData, tvals = t.periods.RDs, intervened_TRT = "gTI.dhigh",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        estimator = "speedglm__glm",
                        family = "quasibinomial")
    tmle_est_dhigh[["estimates"]]

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData,
                        NPMSM = list(surv_dlow, surv_dhigh),
                        wts_data = list(wts.St.dlow, wts.St.dhigh),
                        MSM = MSM.IPAW,
                        AddFUPtables = TRUE,
                        # openFile = FALSE,
                        openFile = TRUE,
                        WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                        file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")
    }

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    MSM.IPAW.both <- survMSM(OData,
                        wts_data = list(dlow = wts.St.dlow, dhigh = wts.St.dhigh),
                        tbreaks = c(1:8,12,16)-1,
                        est_name = "IPAW",
                        getSEs = TRUE,
                        glm_package = "speedglm")

    names(MSM.IPAW.both)
    lapply(MSM.IPAW.both, '[[', 'estimates')

    MSM.IPAW_dlow <- survMSM(OData,
                             wts_data = wts.St.dlow,
                             tbreaks = c(1:8,12,16)-1,
                             est_name = "IPAW",
                             getSEs = TRUE,
                             glm_package = "speedglm")
    names(MSM.IPAW_dlow)

    MSM.IPAW_dhigh <- survMSM(OData,
                              wts_data = wts.St.dhigh,
                              tbreaks = c(1:8,12,16)-1,
                              est_name = "IPAW",
                              getSEs = TRUE,
                              glm_package = "speedglm")


    MSM.IPAW.both[["estimates"]]
    MSM.IPAW_dlow[["estimates"]]
    MSM.IPAW_dhigh[["estimates"]]

    t.periods.RDs <- 0:10
    MSM.RDtables = get_MSM_RDs(MSM.IPAW.both, t.periods.RDs, getSEs = TRUE)
    MSM.RDtables2 = get_RDs(MSM.IPAW.both, "St.MSM", getSEs = TRUE)
    MSM.RDtables
    MSM.RDtables2
}

# ---------------------------------------------------------------------------
# Test xgboost learners
# ---------------------------------------------------------------------------
test.xgboost.10Kdata <- function() {
  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  if (reqxgb) {
    `%+%` <- function(a, b) paste0(a, b)
    # library("stremr")
    options(stremr.verbose = TRUE)
    # options(stremr.verbose = FALSE)
    options(gridisl.verbose = TRUE)
    # options(gridisl.verbose = FALSE)
    require("data.table")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 1,000 IDs
    # Odat_DT <- Odat_DT[ID %in% (1:1000), ]

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
    # FIT PROPENSITY SCORES WITH xgboost glm and no CV
    # ----------------------------------------------------------------
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            estimator = "xgboost__gbm", fit_method = "none",
                            family = "quasibinomial", rounds = 1000, early_stopping_rounds = 50)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    # openFile = FALSE,
                    openFile = TRUE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # model_fits_gA <- OData$modelfit.gA$get.fits()
    # models <- model_fits_gA[[3]]$get_best_models()
    # for (model in models) {
    #     gridisl::print_tables(model)
    # }

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost glm and V fold CV
    # ----------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit_method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(estimator = "xgboost_glm")

    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                            family = "quasibinomial", rounds = 1000, early_stopping_rounds = 50)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # str(OData)
    # library("gridisl")
    # OData$modelfit.gA$get.fits()[[1]]$get_best_models()[[1]]
    # class(OData$modelfit.gA$get.fits()[[1]]$get_best_models()[[1]])

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    # openFile = FALSE,
                    openFile = TRUE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    # openFile = FALSE,
                    openFile = TRUE,
                    ymin = 0.8,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost gbm and no CV
    # ----------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm", fit_method = "none")
    # set_all_stremr_options(estimator = "xgboost_gbm")
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            estimator = "xgboost__gbm", fit_method = "none",
                            family = "quasibinomial", early_stopping_rounds = 10)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)
    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    # openFile = FALSE,
                    openFile = TRUE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost gbm and V fold CV
    # ----------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm", fit_method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(estimator = "xgboost_gbm")

    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                           stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                           estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                           family = "quasibinomial", early_stopping_rounds = 10)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # surv_dlow[["estimates"]][["St.KM"]] <- 1
    pl1 <- ggsurv(list(surv_dlow[["estimates"]], surv_dhigh[["estimates"]]))
    pl2 <- ggsurv(list(surv_dlow[["estimates"]], surv_dhigh[["estimates"]]), surv_name = "St."%+%"KM")

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    openFile = TRUE,
                    # openFile = FALSE,
                    plotKM = TRUE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    openFile = TRUE,
                    # openFile = FALSE,
                    plotKM = TRUE,
                    use_ggplot = FALSE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # OData$Qlearn.fit$get.fits()
    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost glm and CV
    # ---------------------------------------------------------------------------------------------------------
    # params = list(fit.package = "xgboost", fit.algorithm = "glm", family = "quasibinomial") # , objective = "reg:logistic"
    checkException(tmle_est <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE))

    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit_method = "cv", fold_column = "fold_ID")
    t.surv <- c(0:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    tmle_est_dlow <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        estimator = "xgboost__glm",
                        family = "quasibinomial", early_stopping_rounds = 10)
    tmle_est_dlow[["estimates"]]

    tmle_est_dhigh <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        estimator = "xgboost__glm",
                        family = "quasibinomial", early_stopping_rounds = 10)
    tmle_est_dhigh[["estimates"]]

    gcomp_est_dlow <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             estimator = "xgboost__glm",
                             family = "quasibinomial", early_stopping_rounds = 10)
    gcomp_est_dlow[["estimates"]]

    gcomp_est_dhigh <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             estimator = "xgboost__glm",
                             family = "quasibinomial", early_stopping_rounds = 10)
    gcomp_est_dhigh[["estimates"]]

    make_report_rmd(OData,
                    # AddFUPtables = FALSE,
                    openFile = TRUE,
                    # openFile = FALSE,
                    NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    GCOMP = list(gcomp_est_dlow, gcomp_est_dhigh),
                    TMLE = list(tmle_est_dlow, tmle_est_dhigh)
                    )

   # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm", fit_method = "cv", fold_column = "fold_ID")
    t.surv <- c(0:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    tmle_est_dlow <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                        family = "quasibinomial", nrounds  = 10, early_stopping_rounds = 3)
    tmle_est_dlow[["estimates"]]

    tmle_est_dhigh <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                        family = "quasibinomial", nrounds = 10, early_stopping_rounds = 3)
    tmle_est_dhigh[["estimates"]]

    gcomp_est_dlow <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                             family = "quasibinomial", nrounds = 10, early_stopping_rounds = 3)
    gcomp_est_dlow[["estimates"]]

    gcomp_est_dhigh <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                             family = "quasibinomial", nrounds = 10, early_stopping_rounds = 3)
    gcomp_est_dhigh[["estimates"]]

    make_report_rmd(OData,
                # AddFUPtables = FALSE,
                openFile = TRUE,
                # openFile = FALSE,
                NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                GCOMP = list(gcomp_est_dlow, gcomp_est_dhigh),
                TMLE = list(tmle_est_dlow, tmle_est_dhigh)
                )

    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost gbm and CV (WITH EXPLICIT PARAMETER SPECS FOR GBM)
    # ---------------------------------------------------------------------------------------------------------
    # params = list(fit.package = "xgboost", fit.algorithm = "gbm", family = "quasibinomial") # , objective = "reg:logistic"

    t.surv <- c(0:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    params <- gridisl::defModel(estimator = "xgboost__glm", family = "quasibinomial",
                                nthread = 2,
                                nrounds = 10,
                                early_stopping_rounds = 3)

    tmle_est_dlow <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        models = params, fit_method = "cv", fold_column = "fold_ID")

    tmle_est_dhigh <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        models = params, fit_method = "cv", fold_column = "fold_ID")
    tmle_est_dlow[["estimates"]]
    tmle_est_dhigh[["estimates"]]

    gcomp_est_dlow <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             models = params, fit_method = "cv", fold_column = "fold_ID")
    gcomp_est_dhigh <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             models = params, fit_method = "cv", fold_column = "fold_ID")
    gcomp_est_dlow[["estimates"]]
    gcomp_est_dhigh[["estimates"]]

    make_report_rmd(OData,
                # AddFUPtables = FALSE,
                openFile = TRUE,
                # openFile = FALSE,
                NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                GCOMP = list(gcomp_est_dlow, gcomp_est_dhigh),
                TMLE = list(tmle_est_dlow, tmle_est_dhigh)
                )


    # ---------------------------------------------------------------------------------------------------------
    # GCOMP w/ xgboost gbm and CV (WITH EXPLICIT PARAMETER SPECS FOR GBM)
    # ---------------------------------------------------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm", fit_method = "cv", fold_column = "fold_ID")
    # t.surv <- c(0:10)
    t.surv <- c(3)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    params <- gridisl::defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
                                seed = 23,
                                ## learning_rate = 0.1, # learning_rate = 0.01, # learning_rate = 0.05,
                                # nthread = 1,
                                nrounds = 500, early_stopping_rounds = 10,
                                param_grid = list(
                                    learning_rate = c(.1, .3, .5), # .05,
                                    max_depth = c(seq(3, 19, 4), 25),
                                    min_child_weight = c(1, 3, 5, 7),
                                    gamma = c(.0, .05, seq(.1, .9, by=.2), 1),
                                    colsample_bytree = c(.4, .6, .8, 1),
                                    subsample = c(.5, .75, 1),
                                    lambda = c(.1, .5, 2, 5), # lambda = c(1,2,5),
                                    alpha = c(0, .1, .5)
                                    ## Maximum delta step we allow each tree’s weight estimation to be.
                                    ## If the value is set to 0, it means there is no constraint.
                                    ## If it is set to a positive value, it can help making the update step more conservative.
                                    ## Might help in logistic regression when class is extremely imbalanced.
                                    ## Set it to value of 1-10 to help control the update
                                    # max_delta_step = c(0, 1, 2, 5, 10)
                                    )
                                )

    gcomp_est_dlow <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             models = params, fit_method = "cv", fold_column = "fold_ID")

    #   est_name time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type rule.name
    # 1    GCOMP    3 0.9210713      NA          NA          FALSE              4 pooled  gTI.dlow
    # [1] "Surv est: 0.920427939757705"
    # [1] "Surv est 2: 0.920427939757705"
    # [1] "No. of obs for last prediction of Q: 10000"
    # [1] "EY^* estimate at t=3: 0.07957"

    gcomp_est_dhigh <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             models = params, fit_method = "cv", fold_column = "fold_ID")

    gcomp_est_dlow[["estimates"]]
    gcomp_est_dhigh[["estimates"]]

#   est_name time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type rule.name
# 1    GCOMP    3 0.9210713      NA          NA          FALSE              4 pooled gTI.dhigh
# > gcomp_est_dlow[["estimates"]]
#   est_name time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type rule.name
# 1    GCOMP   10 0.7403922      NA          NA          FALSE             11 pooled  gTI.dlow
# > gcomp_est_dhigh[["estimates"]]
#   est_name time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type rule.name
# 1    GCOMP   10 0.7403922      NA          NA          FALSE             11 pooled gTI.dhigh

    make_report_rmd(OData,
                # AddFUPtables = FALSE,
                openFile = TRUE,
                # openFile = FALSE,
                NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
                GCOMP = list(gcomp_est_dlow, gcomp_est_dhigh),
                TMLE = list(tmle_est_dlow, tmle_est_dhigh)
                )

    # ---------------------------------------------------------------------------------------------------------
    # Parallel TMLE w/ xgboost gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------
    # IMPORT DATA
    # ----------------------------------------------------------------
    `%+%` <- function(a, b) paste0(a, b)
    library("stremr")
    library("data.table")
    data.table::setDTthreads(1)
    # library("gridisl")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 1,000 IDs
    # Odat_DT <- Odat_DT[ID %in% (1:1000), ]
    ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
    lagnodes <- c("C", "TI", "N")
    newVarnames <- lagnodes %+% ".tminus1"
    Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
    # indicator that the person has never been on treatment up to current t
    Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
    Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

    setkeyv(Odat_DT, cols = c("ID", "t"))

    OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
    OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    require("foreach")
    require("doParallel")
    stopImplicitCluster()
    registerDoParallel(cores = 4)
    cl <- makeForkCluster(4, outfile = "")
    # cl <- makeForkCluster(10)
    # registerDoParallel(cl)
    # cl <- makeCluster(10)
    # registerDoParallel(cl)
    # data.table::setDTthreads(1)
    # stopCluster(cl)


    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit_method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit_method = "none")

    # params = list(fit.package = "xgboost", fit.algorithm = "gbm", family = "quasibinomial") # , objective = "reg:logistic"
    # params <- gridisl::defModel(estimator = "xgboost__gbm", family = "quasibinomial", nthread = 1,
    #                               nrounds = 500,
    #                               learning_rate = 0.05, # learning_rate = 0.01,
    #                               max_depth = 5,
    #                               min_child_weight = 10,
    #                               colsample_bytree = 0.3,
    #                               early_stopping_rounds = 10,
    #                               seed = 23)

    tmle.model <- "xgb.glm"
    params <- gridisl::defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                # nthread = 2,
                                nrounds = 100,
                                early_stopping_rounds = 5)

    t.surv <- c(1:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    tmle_est <- fitGCOMP(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE, fit_method = "cv", fold_column = "fold_ID",
                        # parallel = FALSE)
                        parallel = TRUE)

    # tmle_est <- fitTMLE(OData, tvals = t.surv,
    #                     intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
    #                     stratifyQ_by_rule = FALSE,
    #                     parallel = TRUE)

    print(tmle_est[["estimates"]])

    # est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates   type
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

    # latest server run using PSOCK cluster:
    #         est_name  t       risk      surv ALLsuccessTMLE nFailedUpdates   type     TMLE_Var     TMLE_SE rule.name
    #  1:     TMLE  1 0.02450699 0.9754930           TRUE              0 pooled 3.674968e-05 0.006062151 gTI.dhigh
    #  2:     TMLE  2 0.04892953 0.9510705           TRUE              0 pooled 7.700189e-05 0.008775072 gTI.dhigh
    #  3:     TMLE  3 0.08690177 0.9130982           TRUE              0 pooled 1.482959e-04 0.012177679 gTI.dhigh
    #  4:     TMLE  4 0.12705143 0.8729486           TRUE              0 pooled 2.157107e-04 0.014687095 gTI.dhigh
    #  5:     TMLE  5 0.14872879 0.8512712           TRUE              0 pooled 2.524471e-04 0.015888583 gTI.dhigh
    #  6:     TMLE  6 0.17629191 0.8237081           TRUE              0 pooled 3.123623e-04 0.017673774 gTI.dhigh
    #  7:     TMLE  7 0.20630433 0.7936957           TRUE              0 pooled 3.562768e-04 0.018875297 gTI.dhigh
    #  8:     TMLE  8 0.22884185 0.7711581           TRUE              0 pooled 3.862008e-04 0.019651993 gTI.dhigh
    #  9:     TMLE  9 0.24857200 0.7514280           TRUE              0 pooled 4.131439e-04 0.020325941 gTI.dhigh
    # 10:     TMLE 10 0.26547898 0.7345210           TRUE              0 pooled 4.322688e-04 0.020791075 gTI.dhigh

    stopCluster(cl)
  }
}


test.xgboost.RFs.10Kdata <- function() {
   # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH randomForest (not implemented)
    # ----------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "drf")
    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                         stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

    # wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    # surv_dlow <- survNPMSM(wts.St.dlow, OData)
    # wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    # surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 openFile = FALSE,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")
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

    # ---------------------------------------------------------------------------------------------------------
    # Define a random grid with xgboost GBMs
    # ---------------------------------------------------------------------------------------------------------
    set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit_method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(estimator = "xgboost_glm")

    model_Grid <- gridisl::defModel(estimator = "xgboost__gbm", family = "quasibinomial",
                                    search_criteria = list(strategy = "RandomDiscrete", max_models = 10),
                                    nrounds = 500, early_stopping_rounds = 10, seed = 123456,
                                    param_grid = list(
                                        eta = c(0.01, 0.02, 0.05, 0.1, 0.3, 0.5),
                                        max_depth = c(4,6,8,10),
                                        max_delta_step = c(0,1),
                                        subsample = c(0.3, 0.5, 0.8, 0.9, 1),
                                        # colsample_bytree = c(0.3, 0.4, 0.5, 0.7, 0.9, 1.0),
                                        colsample_bytree = c(0.5, 0.7, 0.9, 1.0),
                                        scale_pos_weight = 1)
                                    )

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost GBM Grid and V fold CV
    # ----------------------------------------------------------------
    # set_all_stremr_options(estimator = "xgboost_glm")
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            models_CENS = model_Grid, models_TRT = model_Grid, models_MONITOR = model_Grid)

    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                         stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    surv_dlow
    #     t sum_Y_IPAW sum_all_IPAW          ht   St.NPMSM       ht.KM     St.KM rule.name
    # 1   0  1.4491222    159.54285 0.009082965 0.9909170 0.041176471 0.9588235  gTI.dlow
    # 2   1  0.0000000    151.61459 0.000000000 0.9909170 0.000000000 0.9588235  gTI.dlow
    # 3   2  8.1574199    151.64520 0.053792798 0.9376128 0.049079755 0.9117647  gTI.dlow
    # 4   3  1.9198489    136.47299 0.014067611 0.9244229 0.019354839 0.8941176  gTI.dlow
    # 5   4  0.1852475    131.97552 0.001403650 0.9231253 0.006578947 0.8882353  gTI.dlow
    # 6   5  1.7542018    130.94967 0.013396000 0.9107591 0.039735099 0.8529412  gTI.dlow
    # 7   6  0.3535754    124.08692 0.002849417 0.9081640 0.013793103 0.8411765  gTI.dlow
    # 8   7  3.1053348    122.05131 0.025442862 0.8850577 0.027972028 0.8176471  gTI.dlow
    # 9   8  0.0000000    115.64216 0.000000000 0.8850577 0.000000000 0.8176471  gTI.dlow
    # 10  9  2.4588162    115.66551 0.021257990 0.8662431 0.014388489 0.8058824  gTI.dlow
    # 11 10  1.6992984    111.60035 0.015226641 0.8530532 0.021897810 0.7882353  gTI.dlow
    # 12 11  1.4988818    107.51616 0.013940991 0.8411608 0.014925373 0.7764706  gTI.dlow
    # 13 12  3.2758827    104.45601 0.031361361 0.8147808 0.045454545 0.7411765  gTI.dlow
    # 14 13  2.4448476     96.60053 0.025308841 0.7941597 0.023809524 0.7235294  gTI.dlow
    # 15 14  3.4003677     91.93244 0.036987680 0.7647855 0.056910569 0.6823529  gTI.dlow
    # 16 15  2.1977821     83.57311 0.026297717 0.7446734 0.025862069 0.6647059  gTI.dlow
    # 17 16  0.0000000      0.00000         NaN       NaN          NA        NA  gTI.dlow

    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)
    surv_dhigh
    # t sum_Y_IPAW sum_all_IPAW          ht   St.NPMSM       ht.KM     St.KM rule.name
    # 1   0   8.281495     866.8629 0.009553407 0.9904466 0.009324009 0.9906760 gTI.dhigh
    # 2   1  13.018836     780.9117 0.016671328 0.9739345 0.012836970 0.9779587 gTI.dhigh
    # 3   2  18.162327     690.0924 0.026318688 0.9483019 0.023154848 0.9553142 gTI.dhigh
    # 4   3  24.483507     579.4826 0.042250637 0.9082355 0.033388982 0.9234173 gTI.dhigh
    # 5   4  21.954292     477.3127 0.045995618 0.8664606 0.040776699 0.8857634 gTI.dhigh
    # 6   5  10.641143     393.8076 0.027021172 0.8430479 0.024663677 0.8639172 gTI.dhigh
    # 7   6  12.365791     354.5319 0.034879204 0.8136430 0.029702970 0.8382563 gTI.dhigh
    # 8   7  11.293056     307.4146 0.036735592 0.7837534 0.035326087 0.8086439 gTI.dhigh
    # 9   8   7.580146     274.0973 0.027654945 0.7620787 0.032163743 0.7826349 gTI.dhigh
    # 10  9   6.162907     252.5460 0.024403108 0.7434816 0.024844720 0.7631906 gTI.dhigh
    # 11 10   5.236078     235.8163 0.022204054 0.7269733 0.025806452 0.7434953 gTI.dhigh
    # 12 11   4.926659     221.7344 0.022218742 0.7108209 0.020066890 0.7285757 gTI.dhigh
    # 13 12   2.744865     212.1909 0.012935833 0.7016258 0.013698630 0.7185952 gTI.dhigh
    # 14 13   3.181401     205.3523 0.015492409 0.6907559 0.017421603 0.7060761 gTI.dhigh
    # 15 14   3.360858     199.1334 0.016877419 0.6790978 0.017730496 0.6935571 gTI.dhigh
    # 16 15   3.329957     192.4844 0.017299880 0.6673495 0.018050542 0.6810380 gTI.dhigh
    # 17 16   0.000000       0.0000         NaN       NaN          NA        NA gTI.dhigh

    t.surv <- c(10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    gcomp_est <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = model_Grid, stratifyQ_by_rule = FALSE)
    gcomp_est$estimates[]
    # for larger grid (10 random models):
    #    est_name     t      risk      surv ALLsuccessTMLE nFailedUpdates   type rule.name
    #      <char> <num>     <num>     <num>         <lgcl>          <int> <fctr>    <char>
    # 1:    GCOMP    10 0.3287933 0.6712067          FALSE             11 pooled gTI.dhigh

    # for larger grid (20 random models):
    #          <char> <num>     <num>     <num>         <lgcl>          <int> <fctr>    <char>
    # 1:    GCOMP    10 0.3515949 0.6484051          FALSE             11 pooled gTI.dhigh

    tmle_est <- fitTMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = model_Grid, stratifyQ_by_rule = FALSE)
    tmle_est$estimates[]
    ##  RESULTS BASED ON THE LATEST RUN:
    #        est_name     t      risk      surv ALLsuccessTMLE nFailedUpdates   type     TMLE_Var   TMLE_SE rule.name
    #      <char> <num>     <num>     <num>         <lgcl>          <int> <fctr>        <num>     <num>    <char>
    # 1:     TMLE    10 0.2687779 0.7312221           TRUE              0 pooled 0.0004907952 0.0221539 gTI.dhigh

    ## RESULTS BASED ON THE EARLIEST RUNS:
    #    est_name     t       risk      surv ALLsuccessTMLE nFailedUpdates   type TMLE_Var  TMLE_SE rule.name
    #      <char> <num>      <num>     <num>         <lgcl>          <int> <fctr>    <num>    <num>    <char>
    # 1:     TMLE    10 0.05082065 0.9491794           TRUE              0 pooled  6891063 2625.084 gTI.dhigh

    # ---------------------------------------------------------------------------------------------------------
    # Parallel GCOMP w/ xgboost grid search gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    require("doParallel")
    registerDoParallel(cores = 6)
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "glm", fit_method = "cv", fold_column = "fold_ID")
    t.surv <- c(5,6,7,8,9,10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    gcomp_est <- fitGCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE, parallel = FALSE)
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



