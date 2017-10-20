library("stremr")
library("data.table")
setDTthreads(4)

library("foreach")
library("doParallel")
library("xgboost")
library("magrittr")
library("ggplot2")
library("tibble")
library("tidyr")
library("purrr")
library("dplyr")

run_test_xgb_Models <- function(seed){
  data(agaricus.train, package='xgboost')
  data(agaricus.test, package='xgboost')
  dtrain <- xgb.DMatrix(agaricus.train$data, label = agaricus.train$label)
  dtest <- xgb.DMatrix(agaricus.test$data, label = agaricus.test$label)
  watchlist <- list(eval = dtest, train = dtrain)
  param <- list(max_depth = 5, eta = 0.02, nthread = 2, silent = 1,
                objective = "binary:logistic", eval_metric = "auc")
  bst <- xgb.train(param, dtrain, nrounds = 500, watchlist, verbose = FALSE)
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
  fold_column <- "fold_ID"
  fit_method_g <- "cv"

  # ----------------------------------------------------------------
  # FIT PROPENSITY SCORES WITH xgboost gbm and V fold CV
  # ----------------------------------------------------------------
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                          estimator = "xgboost__gbm", fit_method = "cv", fold_column = "fold_ID",
                          family = "quasibinomial", rounds = 1000, early_stopping_rounds = 50)

  # ---------------------------------------------------------------------------------------------------------
  # Parallel test run of xgboost inside function
  # ---------------------------------------------------------------------------------------------------------
  # unregister <- function() {
  #     env <- foreach:::.foreachGlobals
  #     rm(list=ls(name=env), pos=env)
  # }
  # unregister()
  # stopImplicitCluster()

  # cl <- makeForkCluster(4, outfile = "")
  # registerDoParallel(cl); Sys.sleep(2)

  # cat("...running inside run_test_xgb_Models...", "\n")
  # r <- foreach(n=seq.int(8), .packages=c('xgboost'), .export = "run_test_xgb_Models") %dopar% {
  #     run_test_xgb_Models(n)
  # }
  # cat("...finished inside run_test_xgb_Models...", "\n")

  trunc_IPW <- 10
  # tvals <- 0:8
  tvals <- 0:4
  tmax <- 13
  nfolds <- 10 ## number of folds for CV
  # tbreaks = c(1:8,12,16)-1
  tbreaks = c(1:8,11,14)-1
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(tvals)+1))

  ## ------------------------------------------------------------
  ## **** As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** This dataset is to be saved and will be later merged in with all analysis
  ## ------------------------------------------------------------
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh", "gTI.dlow"),
                  trunc_wt = c(FALSE, TRUE),
                  stratifyQ_by_rule = c(TRUE, FALSE)) %>%
                  cross_df() %>%
                  arrange(stratifyQ_by_rule) %>%
                  mutate(nfolds = as.integer(nfolds)) %>%
                  mutate(trunc_MSM = map_dbl(trunc_wt, ~ ifelse(.x, trunc_IPW, Inf))) %>%
                  mutate(trunc_TMLE = trunc_MSM*10)

  ## ------------------------------------------------------------
  ## IPW ANALYSIS
  ## **** For each individual analysis do filter()/subset()/etc to create a grid of parameters specific to given estimator
  ## ------------------------------------------------------------
  IPW_time <- system.time({
    IPW <-  analysis %>%
          rename(trunc_weight = trunc_MSM) %>%
          distinct(intervened_TRT, trunc_weight) %>%

          group_by(intervened_TRT) %>%
          mutate(wts_data = map(first(intervened_TRT), getIPWeights, OData = OData, tmax = tmax)) %>%
          ## save the tables of weights summaries (sep for each regimen)
          mutate(wts_tabs = map(wts_data,
              ~ get_wtsummary(.x, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE))) %>%
          ## save the tables with number at risk / following each rule (sep for each regimen)
          mutate(FUPtimes_tabs = map(wts_data,
                ~ get_FUPtimes(.x, IDnode = ID, tnode = t))) %>%
          ungroup() %>%

          ## IPW-Adjusted KM (Non-Parametric or Saturated MSM):
          mutate(NPMSM = map2(wts_data, trunc_weight,
            ~ survNPMSM(wts_data = .x,
                        trunc_weights = .y,
                        OData = OData))) %>%
          mutate(NPMSM = map(NPMSM, "estimates")) %>%

          ## Crude MSM for hazard (w/out IPW):
          mutate(MSM.crude = map(wts_data,
            ~ survMSM(wts_data = .x,
                      OData = OData,
                      tbreaks = tbreaks,
                      use_weights = FALSE,
                      glm_package = "speedglm"))) %>%
          mutate(MSM.crude = map(MSM.crude, "estimates")) %>%

          ## IPW-MSM for hazard (smoothing over time-intervals in tbreaks):
          mutate(MSM = map2(wts_data, trunc_weight,
            ~ survMSM(wts_data = .x,
                      trunc_weights = .y,
                      OData = OData,
                      tbreaks = tbreaks,
                      glm_package = "speedglm"))) %>%
          mutate(MSM = map(MSM, "estimates")) %>%
          rename(trunc_MSM = trunc_weight)
  })
  IPW_time_hrs <- IPW_time[3]/60/60

  ## save IPW tables (will be later merged with main results dataset)
  IPWtabs <-  analysis %>%
               left_join(IPW) %>%
               distinct(intervened_TRT, trunc_MSM, wts_tabs, FUPtimes_tabs) %>%
               nest(intervened_TRT, wts_tabs, FUPtimes_tabs, .key = "IPWtabs")

  IPW <- IPW %>% select(-wts_data, -wts_tabs, -FUPtimes_tabs)

  ## ------------------------------------------------------------
  ## Parallel GCOMP / TMLE with xgboost gbm and CV
  ## ------------------------------------------------------------
  tmle.model <- "xgb.glm"
  fit_method_Q <- "cv"
  models_Q <- gridisl::defModel(estimator = "xgboost__gbm",
                                family = "quasibinomial",
                                nthread = 2,
                                nrounds = 100,
                                early_stopping_rounds = 20)

  GCOMP_time <- system.time({
    GCOMP <-analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule) %>%
          mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
            ~ fit_GCOMP(intervened_TRT = .x,
                          stratifyQ_by_rule = .y,
                          tvals = tvals,
                          OData = OData,
                          models = models_Q,
                          Qforms = Qforms,
                          fit_method = fit_method_Q,
                          fold_column = fold_column,
                          parallel = TRUE))) %>%
          mutate(GCOMP = map(GCOMP, "estimates"))
  })
  GCOMP_time_hrs <- GCOMP_time[3]/60/60

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  TMLE <- analysis %>%
          rename(trunc_weight = trunc_TMLE) %>%
          distinct(intervened_TRT, stratifyQ_by_rule, trunc_weight)

  TMLE_time <- system.time({
    TMLE <- TMLE %>%
          mutate(TMLE = pmap(TMLE, fit_TMLE,
                             tvals = tvals,
                             OData = OData,
                             models = models_Q,
                             Qforms = Qforms,
                             fit_method = fit_method_Q,
                             fold_column = fold_column,
                             parallel = TRUE)) %>%
          mutate(TMLE = map(TMLE, "estimates")) %>%
          rename(trunc_TMLE = trunc_weight)
  })
  TMLE_time_hrs <- TMLE_time[3]/60/60

  ## ------------------------------------------------------------
  ## COMBINE ALL ANALYSES INTO A SINGLE DATASET
  ## ------------------------------------------------------------
  results <-  analysis %>%
              left_join(IPW) %>%
              left_join(GCOMP) %>%
              left_join(TMLE)

  ## Nest each estimator by treatment regimen (we now only show the main analysis rows)
  results <- results %>%
              # nest(intervened_TRT, NPMSM, MSM.crude, MSM, .key = "estimates")
              nest(intervened_TRT, NPMSM, MSM.crude, MSM, GCOMP, TMLE, .key = "estimates")

  ## Calculate RDs (contrasting all interventions, for each analysis row & estimator).
  ## The RDs data no longer needs the intervened_TRT column
  results <-  results %>%
              mutate(RDs =
                map(estimates,
                  ~ select(.x, -intervened_TRT) %>%
                  map(~ get_RDs(.x)) %>%
                  as_tibble()
                  ))

  ## Clean up by removing the subject-level IC estimates for EVERY SINGLE ESTIMATE / ANALYSIS
  ## WARNING: THIS IS A SIDE-EFFECT FUNCTION!
  # res <- results[["estimates"]] %>%
  #     map(
  #       ~ select(.x, -intervened_TRT) %>%
  #       map(
  #         ~ map(.x,
  #           ~ suppressWarnings(.x[, ("IC.St") := NULL]))))
  # rm(res)

  ## equivalent RD function above but with explicitely defined inside function::
  # results_2 <- results %>%
  #     mutate(RDs = map(estimates, function(.df) {
  #         res <- select(.df, -intervened_TRT) %>%
  #         map(~ get_RDs(.x))
  #         browser()
  #         as_tibble(res)
  #         return(res)
  #         })
  #     )

  ## ------------------------------------------------------------
  ## Add models used for g and Q
  ## Add IPWtabs
  ## ------------------------------------------------------------
  cat("IPW time, hrs: ", IPW_time_hrs, "\n")
  cat("GCOMP time, hrs: ", GCOMP_time_hrs, "\n")
  cat("TMLE time, hrs: ", TMLE_time_hrs, "\n")

  results <- results %>%
             left_join(IPWtabs) %>%
             mutate(fit_method_g = fit_method_g) %>%
             mutate(fit_method_Q = fit_method_Q) %>%
             # mutate(models_g = map(fit_method_g, ~ models_g)) %>%
             mutate(models_Q = map(fit_method_Q, ~ models_Q)) %>%
             mutate(run_time = map(trunc_wt,
              ~ tibble(IPW_time_hrs = IPW_time_hrs, GCOMP_time_hrs = GCOMP_time_hrs, TMLE_time_hrs = TMLE_time_hrs)))

  # as.data.table(results)
  # as.data.table(results)[["models_g"]]
  # as.data.table(results)[["run_time"]]

  return(results)
}

# registerDoParallel(cores = 4); Sys.sleep(2)
cl <- makeForkCluster(3, outfile = "")
registerDoParallel(cl); Sys.sleep(2)

cat("...running outside run_test_xgb_Models...", "\n")
r <- foreach(n=seq.int(8), .packages=c('xgboost')) %dopar% {
    run_test_xgb_Models(n)
}
cat("...finished outside with run_test_xgb_Models...", "\n")

results <- test.xgboost.parallel.10Kdata()

stopCluster(cl)
