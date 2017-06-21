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

    # options(stremr.verbose = TRUE)
    options(stremr.verbose = FALSE)
    # options(gridisl.verbose = TRUE)
    options(gridisl.verbose = FALSE)

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

    ## regularlized glm with h2o
    models_g <<- defModel(estimator = "xgboost__glm", family = "binomial",
                                   nthread = 1,
                                    param_grid = list(
                                        alpha = c(0, 0.5, 1)
                                  ))
    # models_g <<- defModel(estimator = "h2o__glm", family = "binomial",
                                    # nlambdas = 5, lambda_search = TRUE,

    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                           stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                           models_CENS = models_g, models_TRT = models_g,
                           models_MONITOR = defModel(estimator = "speedglm__glm", family = "quasibinomial"),
                          fit_method = "cv", fold_column = "fold_ID"
                          )

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    ## weights based on holdout propensity score predictions
    wts.St.dlow_hold <- getIPWeights(OData, intervened_TRT = "gTI.dlow", holdout = TRUE)
    surv_dlow_hold <- survNPMSM(wts.St.dlow_hold, OData)
    wts.St.dhigh_hold <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", holdout = TRUE)
    surv_dhigh_hold <- survNPMSM(wts.St.dhigh_hold, OData)

    # surv_dlow
    # surv_dlow_hold
    # surv_dhigh
    # surv_dhigh_hold

    # ---------------------------------------------------------------------------------------------------------
    # CV TMLE w/ xgboost gbm and cross-validation selection of Q
    # ---------------------------------------------------------------------------------------------------------
    params <- defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 1,
                                nrounds = 5)

    t.surv <- c(0:2)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    CV_tmle_est <- fit_CVTMLE(OData, tvals = t.surv,
                          intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                          stratifyQ_by_rule = FALSE,
                          fit_method = "cv", # fit_method = "none",
                          fold_column = "fold_ID",
                          parallel = FALSE)
                           # parallel = TRUE)
    # CV_tmle_est[["estimates"]]
    #     est_name time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #  1:     TMLE    1       NA 0.9751527          NA           TRUE              0 pooled 0.001979373
    #  2:     TMLE    2       NA 0.9467156          NA           TRUE              0 pooled 0.003072020
    #  3:     TMLE    3       NA 0.9173932          NA           TRUE              0 pooled 0.003856633
    #  4:     TMLE    4       NA 0.8880877          NA           TRUE              0 pooled 0.004493341
    #  5:     TMLE    5       NA 0.8575970          NA           TRUE              0 pooled 0.005146447
    #  6:     TMLE    6       NA 0.8296318          NA           TRUE              0 pooled 0.005704154
    #  7:     TMLE    7       NA 0.8094443          NA           TRUE              0 pooled 0.005984797
    #  8:     TMLE    8       NA 0.7875892          NA           TRUE              0 pooled 0.006395675
    #  9:     TMLE    9       NA 0.7680555          NA           TRUE              0 pooled 0.006687276
    # 10:     TMLE   10       NA 0.7545343          NA           TRUE              0 pooled 0.006843401
    tmle_est <- fit_TMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        fit_method = "cv", # fit_method = "none",
                        fold_column = "fold_ID",
                        parallel = FALSE)
                        # parallel = TRUE)
    # tmle_est[["estimates"]]
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

    CV_tmle_est <- fit_CVTMLE(OData, tvals = t.surv,
                          intervened_TRT = "gTI.dlow", Qforms = Qforms, models = params,
                          stratifyQ_by_rule = FALSE,
                          fit_method = "cv", # fit_method = "none",
                          fold_column = "fold_ID",
                          parallel = FALSE)
                           # parallel = TRUE)
    # CV_tmle_est[["estimates"]]
    #     est_name time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #  1:     TMLE    1       NA 0.9883379          NA           TRUE              0 pooled 0.001826879
    #  2:     TMLE    2       NA 0.9661213          NA           TRUE              0 pooled 0.004778864
    #  3:     TMLE    3       NA 0.9576946          NA           TRUE              0 pooled 0.005366859
    #  4:     TMLE    4       NA 0.9448527          NA           TRUE              0 pooled 0.006186556
    #  5:     TMLE    5       NA 0.9322860          NA           TRUE              0 pooled 0.006786424
    #  6:     TMLE    6       NA 0.9181835          NA           TRUE              0 pooled 0.007479638
    #  7:     TMLE    7       NA 0.9081946          NA           TRUE              0 pooled 0.007890304
    #  8:     TMLE    8       NA 0.8982760          NA           TRUE              0 pooled 0.008269477
    #  9:     TMLE    9       NA 0.8797293          NA           TRUE              0 pooled 0.008987486
    # 10:     TMLE   10       NA 0.8612292          NA           TRUE              0 pooled 0.009661293
    tmle_est <- fit_TMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dlow", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        fit_method = "cv", # fit_method = "none",
                        fold_column = "fold_ID",
                        parallel = FALSE)
                        # parallel = TRUE)
    # tmle_est[["estimates"]]
    #     est_name time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #  1:     TMLE    1       NA 0.9883911          NA           TRUE              0 pooled 0.001819630
    #  2:     TMLE    2       NA 0.9661918          NA           TRUE              0 pooled 0.004766402
    #  3:     TMLE    3       NA 0.9578060          NA           TRUE              0 pooled 0.005349632
    #  4:     TMLE    4       NA 0.9449523          NA           TRUE              0 pooled 0.006171676
    #  5:     TMLE    5       NA 0.9324055          NA           TRUE              0 pooled 0.006769849
    #  6:     TMLE    6       NA 0.9182947          NA           TRUE              0 pooled 0.007465462
    #  7:     TMLE    7       NA 0.9083105          NA           TRUE              0 pooled 0.007876142
    #  8:     TMLE    8       NA 0.8983794          NA           TRUE              0 pooled 0.008257888
    #  9:     TMLE    9       NA 0.8797762          NA           TRUE              0 pooled 0.008979313
    # 10:     TMLE   10       NA 0.8613221          NA           TRUE              0 pooled 0.009649131

}
