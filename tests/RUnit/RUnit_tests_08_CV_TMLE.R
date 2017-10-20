test.CV_TMLE.10Kdata <- function() {
    reqxgb <- requireNamespace("xgboost", quietly = TRUE)
    reqh2o <- requireNamespace("h2o", quietly = TRUE)
    if (!reqxgb || !reqh2o) return(NULL)

    `%+%` <- function(a, b) paste0(a, b)
    library("stremr")
    library("h2o")
    library("xgboost")
    library("data.table")
    setDTthreads(1)
    library("foreach")
    library("doParallel")

    options(stremr.verbose = FALSE)
    options(gridisl.verbose = FALSE)

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 100 IDs
    # Odat_DT <- Odat_DT[ID %in% (1:100), ]
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
    OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost gbm and V fold CV
    # ----------------------------------------------------------------
    ## xgboost gbm
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                            family = "quasibinomial", rounds = 5, nthread = 1)
    ## h2o gbm
    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                        stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
    #                        estimator = "h2o__gbm", distribution = "bernoulli",
    #                        models_MONITOR = defModel(estimator = "speedglm__glm", family = "quasibinomial"),
    #                       fit_method = "cv", fold_column = "fold_ID", ntrees = 5
    #                       )

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    ## weights based on holdout propensity score predictions
    wts.St.dlow_hold <- getIPWeights(OData, intervened_TRT = "gTI.dlow", holdout = TRUE)
    surv_dlow_hold <- survNPMSM(wts.St.dlow_hold, OData)
    wts.St.dhigh_hold <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", holdout = TRUE)
    surv_dhigh_hold <- survNPMSM(wts.St.dhigh_hold, OData)

    surv_dlow[["St.NPMSM"]]
    surv_dlow_hold[["St.NPMSM"]]
    surv_dhigh[["St.NPMSM"]]
    surv_dhigh_hold[["St.NPMSM"]]

    # ---------------------------------------------------------------------------------------------------------
    # CV TMLE w/ xgboost gbm and cross-validation selection of Q
    # ---------------------------------------------------------------------------------------------------------
    params <- defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 1,
                                nrounds = 5)

    t.surv <- c(0:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    CV_tmle_est <- fit_CVTMLE(OData, tvals = t.surv,
                          intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                          stratifyQ_by_rule = FALSE,
                          fit_method = "cv", # fit_method = "none",
                          fold_column = "fold_ID",
                          parallel = FALSE)
    CV_tmle_est[["estimates"]][["St.TMLE"]]
    CV_tmle_est[["estimates"]][["SE.TMLE"]]

    tmle_est <- fit_TMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        fit_method = "cv", # fit_method = "none",
                        fold_column = "fold_ID",
                        parallel = FALSE)
    tmle_est[["estimates"]][["St.TMLE"]]
    tmle_est[["estimates"]][["SE.TMLE"]]

}
