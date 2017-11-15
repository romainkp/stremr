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
    if (!reqxgb) return(NULL)

    `%+%` <- function(a, b) paste0(a, b)
    library("data.table")
    library("h2o")
    library("stremr")
    library("sl3")
    # options(stremr.verbose = TRUE)
    # options(gridisl.verbose = TRUE)
    options(stremr.verbose = FALSE)
    options(gridisl.verbose = FALSE)
    # set_all_stremr_options(estimator = "speedglm__glm")

    t.periods.RDs <- 0:2
    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 100 IDs
    Odat_DT <- Odat_DT[ID %in% (1:100), ]
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

    # params_g <- gridisl::defModel(estimator = "xgboost__glm", family = "quasibinomial", nthread = 1)
    params_g <- Lrnr_glm_fast$new()
    OData <- fitPropensity(OData,
                            gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            models_CENS = params_g, models_TRT = params_g, models_MONITOR = params_g,
                            fit_method = "cv",
                            fold_column = "fold_ID")

    ## ----------------------------------------------------------------------------------------------------
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.periods.RDs)+1))
    # set_all_stremr_options(estimator = "speedglm__glm")
    params_Q <- Lrnr_glm_fast$new()
    tmle_est_dlow <- fit_TMLE(OData, tvals = t.periods.RDs, intervened_TRT = "gTI.dlow",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        models = params_Q)
    tmle_est_dhigh <- fit_TMLE(OData, tvals = t.periods.RDs, intervened_TRT = "gTI.dhigh",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        models = params_Q)
    tmle_est_dlow[["estimates"]]
    tmle_est_dhigh[["estimates"]]

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
        # make_report_rmd(OData,
        #                 NPMSM = list(surv_dlow, surv_dhigh),
        #                 wts_data = list(wts.St.dlow, wts.St.dhigh),
        #                 MSM = MSM.IPAW,
        #                 AddFUPtables = TRUE,
        #                 # openFile = FALSE,
        #                 openFile = TRUE,
        #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
        #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

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

    # t.periods.RDs <- 0:10
    # MSM.RDtables = get_MSM_RDs(MSM.IPAW.both, t.periods.RDs, getSEs = TRUE)
    # MSM.RDtables2 = get_RDs(MSM.IPAW.both, "St.MSM", getSEs = TRUE)
    # MSM.RDtables
    # MSM.RDtables2
}

# ---------------------------------------------------------------------------
# Test xgboost learners
# ---------------------------------------------------------------------------
test.xgboost.10Kdata <- function() {
    reqxgb <- requireNamespace("xgboost", quietly = TRUE)
    if (!reqxgb) return(NULL)

    `%+%` <- function(a, b) paste0(a, b)
    library("data.table")
    library("h2o")
    library("stremr")
    library("sl3")
    # options(stremr.verbose = TRUE)
    # options(gridisl.verbose = TRUE)
    options(stremr.verbose = FALSE)
    options(gridisl.verbose = FALSE)
    # set_all_stremr_options(estimator = "speedglm__glm")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 100 IDs
    Odat_DT <- Odat_DT[ID %in% (1:100), ]
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

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost glm and no CV
    # ----------------------------------------------------------------
    models_TRT <- Lrnr_xgboost$new(objective = "reg:logistic", rounds = 5, nthread = 1)
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            models_TRT = models_TRT,
                            fit_method = "none")

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 openFile = FALSE,
    #                 # openFile = TRUE,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

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

    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                         stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
    #                         estimator = "xgboost__glm", fit_method = "cv", fold_column = "fold_ID",
    #                         family = "quasibinomial", nthread = 1)

    # wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    # surv_dlow <- survNPMSM(wts.St.dlow, OData)
    # wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    # surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # str(OData)
    # library("gridisl")
    # OData$modelfit.gA$get.fits()[[1]]$get_best_models()[[1]]
    # class(OData$modelfit.gA$get.fits()[[1]]$get_best_models()[[1]])

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 openFile = FALSE,
    #                 # openFile = TRUE,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 openFile = FALSE,
    #                 # openFile = TRUE,
    #                 ymin = 0.8,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost gbm and V fold CV
    # ----------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm", fit_method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(estimator = "xgboost_gbm")
    models_TRT <- Lrnr_xgboost$new(objective = "reg:logistic", rounds = 5, nthread = 1)

    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                           stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                           models_TRT = models_TRT,
                           fit_method = "cv", fold_column = "fold_ID")

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # surv_dlow[["estimates"]][["St.KM"]] <- 1
    pl1 <- ggsurv(list(surv_dlow[["estimates"]], surv_dhigh[["estimates"]]))
    pl2 <- ggsurv(list(surv_dlow[["estimates"]], surv_dhigh[["estimates"]]), surv_name = "St."%+%"KM")

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 # openFile = TRUE,
    #                 openFile = FALSE,
    #                 plotKM = TRUE,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    #     make_report_rmd(OData, NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 AddFUPtables = TRUE,
    #                 # openFile = TRUE,
    #                 openFile = FALSE,
    #                 plotKM = TRUE,
    #                 use_ggplot = FALSE,
    #                 WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # OData$Qlearn.fit$get.fits()
    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost glm and CV
    # ---------------------------------------------------------------------------------------------------------
    # params = list(fit.package = "xgboost", fit.algorithm = "glm", family = "quasibinomial") # , objective = "reg:logistic"
    t.surv <- c(0:2)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    params <- Lrnr_xgboost$new(objective = "reg:logistic", rounds = 5, nthread = 1)
    tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    # make_report_rmd(OData,
    #                 # AddFUPtables = FALSE,
    #                 # openFile = TRUE,
    #                 openFile = FALSE,
    #                 NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #                 GCOMP = list(gcomp_est_dlow, gcomp_est_dhigh),
    #                 TMLE = list(tmle_est_dlow, tmle_est_dhigh)
    #                 )

   # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ xgboost gbm and CV
    # ---------------------------------------------------------------------------------------------------------
    # set_all_stremr_options(fit.package = "xgboost", fit.algorithm = "gbm", fit_method = "cv", fold_column = "fold_ID")
    t.surv <- c(0:2)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    params <- Lrnr_xgboost$new(objective = "reg:logistic", rounds = 5, nthread = 1)
    tmle_est_dlow <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                        Qforms = Qforms, stratifyQ_by_rule = FALSE,
                        models = params,
                        fit_method = "cv", fold_column = "fold_ID")
    # tmle_est_dhigh <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
    #                     Qforms = Qforms, stratifyQ_by_rule = FALSE,
    #                     estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
    #                     family = "quasibinomial", nrounds = 5, early_stopping_rounds = 2, nthread = 1)
    tmle_est_dlow[["estimates"]]
    # tmle_est_dhigh[["estimates"]]

    gcomp_est_dlow <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             models = params, 
                             fit_method = "cv", fold_column = "fold_ID")
    # gcomp_est_dhigh <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh",
    #                          Qforms = Qforms, stratifyQ_by_rule = FALSE,
    #                          estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
    #                          family = "quasibinomial", nrounds = 5, early_stopping_rounds = 2, nthread = 1)
    gcomp_est_dlow[["estimates"]]
    # gcomp_est_dhigh[["estimates"]]

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    # make_report_rmd(OData,
    #             # AddFUPtables = FALSE,
    #             # openFile = TRUE,
    #             openFile = FALSE,
    #             NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #             GCOMP = list(gcomp_est_dlow, gcomp_est_dhigh),
    #             TMLE = list(tmle_est_dlow, tmle_est_dhigh)
    #             )

    # ---------------------------------------------------------------------------------------------------------
    # GCOMP w/ xgboost gbm and CV (WITH EXPLICIT PARAMETER SPECS FOR GBM)
    # ---------------------------------------------------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    t.surv <- c(2)
    # t.surv <- c(0:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    params <- gridisl::defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                search_criteria = list(strategy = "RandomDiscrete", max_models = 3),
                                seed = 23,
                                nthread = 1,
                                nrounds = 5, early_stopping_rounds = 2,
                                param_grid = list(
                                    learning_rate = c(.1, .3, .5), # .05,
                                    max_depth = c(seq(3, 19, 4), 25),
                                    min_child_weight = c(1, 3, 5, 7),
                                    gamma = c(.0, .05, seq(.1, .9, by=.2), 1),
                                    colsample_bytree = c(.4, .6, .8, 1),
                                    subsample = c(.5, .75, 1),
                                    lambda = c(.1, .5, 2, 5), # lambda = c(1,2,5),
                                    alpha = c(0, .1, .5)
                                    )
                                )


    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                           stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)
    tmle_est_dlow <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow",
                             Qforms = Qforms, stratifyQ_by_rule = FALSE,
                             models = params, fit_method = "cv", fold_column = "fold_ID")
    tmle_est_dlow[["estimates"]]

    # if (rmarkdown::pandoc_available(version = "1.12.3"))
    # make_report_rmd(OData,
    #             # AddFUPtables = FALSE,
    #             # openFile = TRUE,
    #             openFile = FALSE,
    #             NPMSM = list(surv_dlow, surv_dhigh), wts_data = list(wts.St.dlow, wts.St.dhigh),
    #             TMLE = list(tmle_est_dlow)
    #             )
}
