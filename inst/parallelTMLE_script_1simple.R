library("xgboost")
library("data.table")
setDTthreads(4)
library("foreach")
library("doParallel")
library("gridisl")
library("stremr")
library("xgboost")

run_test_xgb_Models <- function(seed){
    data(agaricus.train, package='xgboost')
    data(agaricus.test, package='xgboost')
    dtrain <- xgb.DMatrix(agaricus.train$data, label = agaricus.train$label)
    dtest <- xgb.DMatrix(agaricus.test$data, label = agaricus.test$label)
    watchlist <- list(eval = dtest, train = dtrain)
    param <- list(max_depth = 5, eta = 0.02, nthread = 1, silent = 0,
                  objective = "binary:logistic", eval_metric = "auc")
    bst <- xgb.train(param, dtrain, nrounds = 500, watchlist)
    return(bst)
}

test.xgboost.parallel.10Kdata <- function() {
    `%+%` <- function(a, b) paste0(a, b)
    # options(stremr.verbose = TRUE)
    options(stremr.verbose = FALSE)
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
    # Parallel test run of xgboost inside function
    # ---------------------------------------------------------------------------------------------------------
    # cl <- makeForkCluster(4, outfile = "")
    # registerDoParallel(cl); Sys.sleep(2)
    cat("...running inside run_test_xgb_Models...", "\n")
    r <- foreach(n=seq.int(8), .packages=c('xgboost'), .export = "run_test_xgb_Models") %dopar% {
        run_test_xgb_Models(n)
    }
    cat("...finished inside run_test_xgb_Models...", "\n")

    # ---------------------------------------------------------------------------------------------------------
    # Parallel TMLE w/ xgboost gbm and CV
    # ---------------------------------------------------------------------------------------------------------

    tmle.model <- "xgb.glm"
    params <- gridisl::defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 1,
                                nrounds = 100,
                                early_stopping_rounds = 20)

    t.surv <- c(1:10)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    tmle_est <- fit_GCOMP(OData, tvals = t.surv,
                        intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
                        stratifyQ_by_rule = FALSE,
                        # fit_method = "none",
                        fit_method = "cv",
                        fold_column = "fold_ID",
                        parallel = FALSE)
                        # parallel = TRUE)


}

# unregister <- function() {
#     env <- foreach:::.foreachGlobals
#     rm(list=ls(name=env), pos=env)
# }
# unregister()
# stopImplicitCluster()

# registerDoParallel(cores = 4); Sys.sleep(2)
cl <- makeForkCluster(4, outfile = "")
registerDoParallel(cl); Sys.sleep(2)

cat("...running outside run_test_xgb_Models...", "\n")
r <- foreach(n=seq.int(8), .packages=c('xgboost')) %dopar% {
    run_test_xgb_Models(n)
}
cat("...finished outside with run_test_xgb_Models...", "\n")

test.xgboost.parallel.10Kdata()

stopCluster(cl)
