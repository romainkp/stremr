test.iTMLE.10Kdata <- function() {
    reqxgb <- requireNamespace("xgboost", quietly = TRUE)
    reqh2o <- requireNamespace("h2o", quietly = TRUE)
    origami <- requireNamespace("origami", quietly = TRUE)
    if (!reqxgb || !reqh2o || !origami) return(NULL)

    library("origami")
    library("h2o")
    library("xgboost")
    library("data.table")
    library("foreach")
    library("doParallel")
    library("stremr")

    if (length(find("origami_SuperLearner"))==0) {
      return(NULL)
    }

    data.table::setDTthreads(1)
    options(gridisl.verbose = FALSE)
    options(stremr.verbose = FALSE)

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 100 IDs
    Odat_DT <- Odat_DT[ID %in% (1:200), ]
    setkeyv(Odat_DT, cols = c("ID", "t"))

    # ---------------------------------------------------------------------------
    # Define some summaries (lags C[t-1], A[t-1], N[t-1])
    # ---------------------------------------------------------------------------
    ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
    lagnodes <- c("C", "TI", "N")
    newVarnames <- paste0(lagnodes, ".tminus1")
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
    params <- defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 1,
                                nrounds = 5,
                                early_stopping_rounds = 2)

    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, 
                            models_TRT = params,
                            models_CENS = params,
                            # gform_MONITOR = gform_MONITOR,
                            # estimator = "xgboost__gbm",
                            fit_method = "cv",
                            fold_column = "fold_ID")

    # ## h2o gbm
    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                        stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
    #                        estimator = "h2o__gbm", distribution = "bernoulli",
    #                        models_MONITOR = defModel(estimator = "speedglm__glm", family = "quasibinomial"),
    #                       fit_method = "cv", fold_column = "fold_ID"
    #                       )
    # ## regularlized glm with h2o
    # models_g <- defModel(estimator = "h2o__glm", family = "binomial",
    #                                 nlambdas = 5, lambda_search = TRUE,
    #                                 param_grid = list(
    #                                     alpha = c(0, 0.5, 1)
    #                               ))
    # OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
    #                        stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
    #                        models_CENS = models_g, models_TRT = models_g,
    #                        models_MONITOR = defModel(estimator = "speedglm__glm", family = "quasibinomial"),
    #                       fit_method = "cv", fold_column = "fold_ID"
    #                       )

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv_dlow <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv_dhigh <- survNPMSM(wts.St.dhigh, OData)

    # ---------------------------------------------------------------------------------------------------------
    # CV TMLE w/ xgboost gbm and cross-validation selection of Q
    # ---------------------------------------------------------------------------------------------------------
    params <- defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 1,
                                nrounds = 5,
                                early_stopping_rounds = 2)

    # t.surv <- c(1:10)
    t.surv <- 4
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

    SDR_est <- fit_iTMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        # fit_method = "cv",
                        fit_method = "none",
                        # fold_column = "fold_ID",
                        parallel = FALSE,
                        return_fW = TRUE)
                        # parallel = TRUE)
    SDR_est[["estimates"]]
    fW_fit <- SDR_est[["estimates"]][["fW_fit"]][[1]]
    preds_fW <- gridisl::predict_SL(fW_fit, Odat_DT[t==0, ])
    MSE_err <- mean((Odat_DT[t==0, ][["Y.tplus1"]] - preds_fW)^2)
    # MSE_err # [1] 0.07085635

    #         est_name  time    St.SDR   type rule.name
    #       <char> <int>     <num> <char>    <char>
    #  1:    SeqDR     1 0.9752629 pooled gTI.dhigh
    #  2:    SeqDR     2 0.9469588 pooled gTI.dhigh
    #  3:    SeqDR     3 0.9173074 pooled gTI.dhigh
    #  4:    SeqDR     4 0.8877367 pooled gTI.dhigh
    #  5:    SeqDR     5 0.8575215 pooled gTI.dhigh
    #  6:    SeqDR     6 0.8309006 pooled gTI.dhigh
    #  7:    SeqDR     7 0.8101947 pooled gTI.dhigh
    #  8:    SeqDR     8 0.7883888 pooled gTI.dhigh
    #  9:    SeqDR     9 0.7699733 pooled gTI.dhigh
    # 10:    SeqDR    10 0.7559333 pooled gTI.dhigh

    ## WILL RETURN ERROR WITH REGULAR VERSION OF XGBOOST
    ## xgboost has to be recompiled to allow logistic loss y>1 or y<0.
    # DR_trans_est <- fit_iTMLE(OData, tvals = t.surv,
    #                        intervened_TRT = "gTI.dhigh", Qforms = Qforms,
    #                        stratifyQ_by_rule = FALSE,
    #                        fit_method = "none",
    #                        models = params,
    #                        return_fW = TRUE,
    #                        use_DR_transform = TRUE # stabilize = FALSE,
    #                       )
    # DR_trans_est[["estimates"]]

    tmle_est <- fit_TMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        fit_method = "cv", # fit_method = "none",
                        fold_column = "fold_ID",
                        parallel = FALSE,
                        return_fW = TRUE)
                        # parallel = TRUE)
    tmle_est[["estimates"]]
    # fW_fit <- tmle_est[["estimates"]][["fW_fit"]][[1]]
    # preds_fW <- gridisl::predict_SL(fW_fit, Odat_DT[t==0, ])
    # MSE_err <- mean((Odat_DT[t==0, ][["Y.tplus1"]] - preds_fW)^2)
    # MSE_err # [1] 0.07085635
    #     est_name  time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #       <char> <int>   <lgcl>     <num>      <lgcl>         <lgcl>          <int> <char>       <num>
    #  1:     TMLE     1       NA 0.9751499          NA           TRUE              0 pooled 0.001968403
    #  2:     TMLE     2       NA 0.9466908          NA           TRUE              0 pooled 0.003061786
    #  3:     TMLE     3       NA 0.9173085          NA           TRUE              0 pooled 0.003840221
    #  4:     TMLE     4       NA 0.8880466          NA           TRUE              0 pooled 0.004475643
    #  5:     TMLE     5       NA 0.8575049          NA           TRUE              0 pooled 0.005122649
    #  6:     TMLE     6       NA 0.8294285          NA           TRUE              0 pooled 0.005668354
    #  7:     TMLE     7       NA 0.8094311          NA           TRUE              0 pooled 0.005952984
    #  8:     TMLE     8       NA 0.7879003          NA           TRUE              0 pooled 0.006331397
    #  9:     TMLE     9       NA 0.7684924          NA           TRUE              0 pooled 0.006624978
    # 10:     TMLE    10       NA 0.7549844          NA           TRUE              0 pooled 0.006795919

    SDR_est <- fit_iTMLE(OData, tvals = t.surv,
                          intervened_TRT = "gTI.dlow", Qforms = Qforms, models = params,
                          stratifyQ_by_rule = FALSE,
                          # fit_method = "cv",
                          fit_method = "none",
                          fold_column = "fold_ID",
                          parallel = FALSE)
                           # parallel = TRUE)
    # SDR_est[["estimates"]]
    #       <char> <int>     <num> <char>    <char>
    #  1:    SeqDR     1 0.9882525 pooled  gTI.dlow
    #  2:    SeqDR     2 0.9663698 pooled  gTI.dlow
    #  3:    SeqDR     3 0.9579345 pooled  gTI.dlow
    #  4:    SeqDR     4 0.9451448 pooled  gTI.dlow
    #  5:    SeqDR     5 0.9327012 pooled  gTI.dlow
    #  6:    SeqDR     6 0.9183899 pooled  gTI.dlow
    #  7:    SeqDR     7 0.9082106 pooled  gTI.dlow
    #  8:    SeqDR     8 0.8985577 pooled  gTI.dlow
    #  9:    SeqDR     9 0.8797818 pooled  gTI.dlow
    # 10:    SeqDR    10 0.8615536 pooled  gTI.dlow

    tmle_est <- fit_TMLE(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dlow", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        fit_method = "cv", # fit_method = "none",
                        fold_column = "fold_ID",
                        parallel = FALSE)
                        # parallel = TRUE)
    # tmle_est[["estimates"]]
    #     est_name  time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
    #       <char> <int>   <lgcl>     <num>      <lgcl>         <lgcl>          <int> <char>       <num>
    #  1:     TMLE     1       NA 0.9883911          NA           TRUE              0 pooled 0.001819630
    #  2:     TMLE     2       NA 0.9661918          NA           TRUE              0 pooled 0.004766402
    #  3:     TMLE     3       NA 0.9578060          NA           TRUE              0 pooled 0.005349632
    #  4:     TMLE     4       NA 0.9449523          NA           TRUE              0 pooled 0.006171676
    #  5:     TMLE     5       NA 0.9324055          NA           TRUE              0 pooled 0.006769849
    #  6:     TMLE     6       NA 0.9182947          NA           TRUE              0 pooled 0.007465462
    #  7:     TMLE     7       NA 0.9083105          NA           TRUE              0 pooled 0.007876142
    #  8:     TMLE     8       NA 0.8983794          NA           TRUE              0 pooled 0.008257888
    #  9:     TMLE     9       NA 0.8797762          NA           TRUE              0 pooled 0.008979313
    # 10:     TMLE    10       NA 0.8613221          NA           TRUE              0 pooled 0.009649131

}
