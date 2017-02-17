
test.CV_TMLE.10Kdata <- function() {
    `%+%` <- function(a, b) paste0(a, b)
    library("xgboost")
    library("data.table")
    setDTthreads(1)
    library("foreach")
    library("doParallel")

    library("gridisl")
    # library("stremr")
    library("xgboost")
    library("data.table")
    setDTthreads(1)

    options(stremr.verbose = TRUE)
    # options(stremr.verbose = FALSE)
    options(gridisl.verbose = TRUE)
    # options(gridisl.verbose = FALSE)

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
    OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    OData <- define_CVfolds(OData, nfolds = 5, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH xgboost gbm and V fold CV
    # ----------------------------------------------------------------
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                            family = "quasibinomial", rounds = 1000, early_stopping_rounds = 50)

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # ---------------------------------------------------------------------------------------------------------
    # CV TMLE w/ xgboost gbm and cross-validation selection of Q
    # ---------------------------------------------------------------------------------------------------------
    tmle.model <- "xgb.glm"
    params <- gridisl::defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 2,
                                nrounds = 100,
                                early_stopping_rounds = 2)

    t.surv <- c(1:10)
    Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    CV_tmle_est <- fitCVTMLE(OData, tvals = t.surv,
                          intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                          stratifyQ_by_rule = FALSE,
                          fit_method = "cv", # fit_method = "none",
                          fold_column = "fold_ID",
                          parallel = FALSE)
                           # parallel = TRUE)
    CV_tmle_est[["estimates"]]
    #     est_name time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #  1:     TMLE    1       NA 0.9751897          NA           TRUE              0 pooled 0.001976476
    #  2:     TMLE    2       NA 0.9467307          NA           TRUE              0 pooled 0.003070247
    #  3:     TMLE    3       NA 0.9174615          NA           TRUE              0 pooled 0.003850820
    #  4:     TMLE    4       NA 0.8881408          NA           TRUE              0 pooled 0.004487804
    #  5:     TMLE    5       NA 0.8576092          NA           TRUE              0 pooled 0.005140892
    #  6:     TMLE    6       NA 0.8296722          NA           TRUE              0 pooled 0.005690897
    #  7:     TMLE    7       NA 0.8094901          NA           TRUE              0 pooled 0.005973558
    #  8:     TMLE    8       NA 0.7877354          NA           TRUE              0 pooled 0.006376092
    #  9:     TMLE    9       NA 0.7681367          NA           TRUE              0 pooled 0.006671797
    # 10:     TMLE   10       NA 0.7546043          NA           TRUE              0 pooled 0.006828314

    tmle_est <- fitTMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        fit_method = "cv", # fit_method = "none",
                        fold_column = "fold_ID",
                        parallel = FALSE)
                        # parallel = TRUE)
    #     tmle_est[["estimates"]]
    #     est_name time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #  1:     TMLE    1       NA 0.9751499          NA           TRUE              0 pooled 0.001968397
    #  2:     TMLE    2       NA 0.9466906          NA           TRUE              0 pooled 0.003061770
    #  3:     TMLE    3       NA 0.9173085          NA           TRUE              0 pooled 0.003840221
    #  4:     TMLE    4       NA 0.8880466          NA           TRUE              0 pooled 0.004475643
    #  5:     TMLE    5       NA 0.8575049          NA           TRUE              0 pooled 0.005122649
    #  6:     TMLE    6       NA 0.8293748          NA           TRUE              0 pooled 0.005665724
    #  7:     TMLE    7       NA 0.8094311          NA           TRUE              0 pooled 0.005952984
    #  8:     TMLE    8       NA 0.7880349          NA           TRUE              0 pooled 0.006332490
    #  9:     TMLE    9       NA 0.7684924          NA           TRUE              0 pooled 0.006624978
    # 10:     TMLE   10       NA 0.7549844          NA           TRUE              0 pooled 0.006795919

}
