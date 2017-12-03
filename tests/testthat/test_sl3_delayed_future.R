context("sl3 with delayed future")

  # delayed_object_7 <- delayed(rnorm(1000000))
  # delayed_object_3 <- delayed(rnorm(1000000))
  # delayed_object_6 <- delayed(rnorm(1000000))
  # delayed_object_7 <- delayed(rnorm(1000000))
  # bundle <- bundle_delayed(list(delayed_object_7, delayed_object_3, delayed_object_6, delayed_object_7))
  # adder <- function(x, y){x + y}
  # delayed_adder <- delayed_fun(adder)
  # chained_delayed_10 <- delayed_adder(delayed_object_7, delayed_object_3)
  # # compute it using the future plan (two multicore workers), verbose mode lets us
  # # see the computation order
  # res <- bundle$compute(nworkers = 2)
  # res <- chained_delayed_10$compute(nworkers = 2)

  ## -----------------------------------------------------------------------
  ## Analyses by intervention
  ## -----------------------------------------------------------------------
  # devtools::install_github("osofr/stremr", ref = "delayed_future", dependencies = FALSE)
  # devtools::install_github("jeremyrcoyle/sl3")
  # devtools::install_github("jeremyrcoyle/delayed")
  # devtools::install_github("jeremyrcoyle/origami")
  # plan(multicore, workers = 4)
  # plan(multisession)
  # plan(sequential)
  # library(delayed)

  library("SuperLearner")
  library("future")
  library("delayed")
  library("sl3")
  library("stremr")
  library("data.table")
  library("magrittr")
  library("ggplot2")
  library("tibble")
  library("tidyr")
  library("purrr")
  library("dplyr")

  data.table::setDTthreads(1)
  
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  options(sl3.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = TRUE)
  # options(sl3.verbose = TRUE)

  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # select only the first 1000 IDs
  Odat_DT <- Odat_DT[ID %in% (1:500), ]
  setkeyv(Odat_DT, cols = c("ID", "t"))

  ## -----------------------------------------------------------------------
  ## Define some summaries (lags C[t-1], A[t-1], N[t-1])
  ## -----------------------------------------------------------------------
  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- paste0(lagnodes, ".tminus1")
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

  ## ------------------------------------------------------------------
  ## Propensity score models for Treatment, Censoring & Monitoring
  ## ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))

  ## ------------------------------------------------------------
  ## **** As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** This dataset is to be saved and will be later merged in with all analysis
  ## ------------------------------------------------------------
  # tvals <- 0:10
  tvals <- 2
  ## This dataset defines all parameters that we like to vary in this analysis (including different interventions)
  ## That is, each row of this dataset corresponds with a single analysis, for one intervention of interest.
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh"),
                  stratifyQ_by_rule = c(FALSE)) %>%
                  cross_df() %>%
                  arrange(stratifyQ_by_rule)

  ## ------------------------------------------------------------------------
  ## Define models for fitting propensity scores (g) -- PARAMETRIC LOGISTIC REGRESSION
  ## ------------------------------------------------------------------------
  fit_method_g <- "none"
  # models_g <- defModel(estimator = "speedglm__glm", family = "quasibinomial")
  lrn_glm <- Lrnr_glm_fast$new()
  lrn_glm_sm <- Lrnr_glm_fast$new(covariates = c("CVD"))
  lrn_glmnet_binom <- Lrnr_pkg_SuperLearner$new("SL.glmnet")
  lrn_glmnet_gaus <- Lrnr_pkg_SuperLearner$new("SL.glmnet")
  sl <- Lrnr_sl$new(learners = Stack$new(lrn_glm, lrn_glm_sm), # , lrn_glmnet_binom
                    metalearner = Lrnr_nnls$new())
  models_g <- sl

  OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), TRT = "TI", OUTCOME = "Y.tplus1") %>%
           stremr::define_CVfolds(nfolds = 10, fold_column = "fold_ID")

  # subset_idx <- which(!is.na(OData$dat.sVar[["TI"]]) & !is.na(OData$dat.sVar[["highA1c"]]) & !is.na(OData$dat.sVar[["lastNat1"]]))

  # fit <- stremr:::test_sl3_fit_single_regression(OData,
  #         OData$nodes,
  #         models = models_g,
  #         predvars = c("highA1c", "lastNat1"),
  #         outvar = "TI",
  #         subset_idx = subset_idx)

  # plan(multicore)
  # fit <- stremr:::test_sl3_fit_single_regression(OData,
  #         OData$nodes,
  #         models = models_g,
  #         predvars = c("highA1c", "lastNat1"),
  #         outvar = "TI",
  #         subset_idx = subset_idx)

  # task <- sl3::sl3_Task$new(Odat_DT[!is.na(TI) & !is.na(highA1c) & !is.na(lastNat1),],
  #                           covariates = c("highA1c", "lastNat1"),
  #                           outcome = "TI",
  #                           id = "ID")
  # tmc <- system.time(model.fit <- models_g$train(task))


  # plan(sequential)
  # task <- sl3::sl3_Task$new(Odat_DT[!is.na(TI) & !is.na(highA1c) & !is.na(lastNat1),],
  #                           covariates = c("highA1c", "lastNat1"),
  #                           outcome = "TI",
  #                           id = "ID")
  # tseq <- system.time(model.fit <- models_g$train(task))

  # plan(multisession)
  # task <- sl3::sl3_Task$new(Odat_DT[!is.na(TI) & !is.na(highA1c) & !is.na(lastNat1),],
  #                           covariates = c("highA1c", "lastNat1"),
  #                           outcome = "TI",
  #                           id = "ID")
  # tmsesh <- system.time(model.fit <- models_g$train(task))
  # print(tmsesh)

  ## ------------------------------------------------------------------------
  ## Define models for iterative G-COMP (Q) -- PARAMETRIC LOGISTIC REGRESSION
  ## ------------------------------------------------------------------------
  ## regression formulas for Q's:
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + lastNat1 + TI + TI.tminus1", (max(tvals)+1))
  ## no cross-validation model selection, just fit a single model specified below
  fit_method_Q <- "none"
  ## Use speedglm to fit all Q.
  ## NOTE: it is currently not possible to use fit_method_Q <- "cv" with speedglm or glm.
  ## To perform cross-validation with GLM use 'estimator="h2o__glm"' or 'estimator="xgboost__glm"'
  # models_Q <- defModel(estimator = "speedglm__glm", family = "quasibinomial")
  # models_Q <- defModel(estimator = "xgboost__glm", family = "quasibinomial")
  models_Q <- defModel(estimator = "xgboost__gbm", family = "quasibinomial", nrounds = 5, nthread = 5)
  # models_Q <- lrn_glm
  # models_Q <- sl

  ## ----------------------------------------------------------------
  ## Fit propensity score models.
  ## We are using the same model ensemble defined in models_g for censoring, treatment and monitoring mechanisms.
  ## ----------------------------------------------------------------

  # plan(sequential)
  # t_run_seq <- system.time(
    OData <- fitPropensity(OData,
                            gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT,
                            models_TRT = models_g,
                            fit_method = fit_method_g
                            )
    # )
  # print("t_run_seq: "); print(t_run_seq)

  # plan(multisession, workers = 5)
  # t_run_multisession <- system.time(
  #   OData <- fitPropensity(OData,
  #                           gform_TRT = gform_TRT,
  #                           stratify_TRT = stratify_TRT,
  #                           models_TRT = models_g,
  #                           fit_method = fit_method_g
  #                           )
  #   )
  # print("t_run_multisession: "); print(t_run_multisession)

  # plan(multicore, workers = 5)
  # t_run_multicore <- system.time(
  #   OData <- fitPropensity(OData,
  #                           gform_TRT = gform_TRT,
  #                           stratify_TRT = stratify_TRT,
  #                           models_TRT = models_g,
  #                           fit_method = fit_method_g
  #                           )
  #   )
  # print("t_run_multicore"); print(t_run_multicore)
## ERROR:
# Assertion failure at kmp_runtime.cpp(6480): __kmp_thread_pool == __null.
# OMP: Error #13: Assertion failure at kmp_runtime.cpp(6480).
# OMP: Hint: Please submit a bug report with this message, compile and run commands used, and machine configuration info including native compiler and operating system versions. Faster response will be obtained by including all program sources. For information on submitting this issue, please see http://www.intel.com/software/products/support/.


  ## Get the dataset with weights:
  wts_data <- getIPWeights(intervened_TRT = "gTI.dlow", OData = OData)

  ## ------------------------------------------------------------
  ## Parallel GCOMP ANALYSIS
  ## ------------------------------------------------------------
  # plan(multicore, workers = 10)
  # plan(multisession, workers = 10)
  GCOMP <-analysis %>%
        distinct(intervened_TRT, stratifyQ_by_rule) %>%
        mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
          ~ fit_GCOMP(intervened_TRT = .x,
                        stratifyQ_by_rule = .y,
                        tvals = tvals,
                        OData = OData,
                        models = models_Q,
                        Qforms = Qforms,
                        fit_method = fit_method_Q
                        ))) %>%
        mutate(GCOMP = map(GCOMP, "estimates"))


  # test_that("GCOMP results w/out Monitoring match", {
  #   GCOMP[["GCOMP"]][[1]][["St.GCOMP"]]
  #   # [1] 0.99 0.99 0.99

  #   GCOMP[["GCOMP"]][[2]][["St.GCOMP"]]
  #   # [1] 0.9900000 0.9719893 0.9593532
  # })

## with rare outcomes
# [[1]]
#    est_name time  St.GCOMP St.TMLE ALLsuccessTMLE nFailedUpdates       type              IC.St fW_fit rule.name
# 1:    GCOMP   10 0.9828276      NA          FALSE             11 stratified NA,NA,NA,NA,NA,NA,   NULL  gTI.dlow
# [[2]]
#    est_name time  St.GCOMP St.TMLE ALLsuccessTMLE nFailedUpdates       type              IC.St fW_fit rule.name
# 1:    GCOMP   10 0.9778367      NA          FALSE             11 stratified NA,NA,NA,NA,NA,NA,   NULL gTI.dhigh

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  TMLE <- CVTMLE <- analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule)

  TMLE <- TMLE %>%
        mutate(TMLE = pmap(TMLE, fit_TMLE,
                           tvals = tvals,
                           OData = OData,
                           models = models_Q,
                           Qforms = Qforms,
                           fit_method = fit_method_Q
                           )) %>%
        mutate(TMLE = map(TMLE, "estimates"))

  # test_that("TMLE results w/out Monitoring match", {
  #   TMLE[["TMLE"]][[1]][["St.TMLE"]]
  #   # [1] 0.99 0.99 0.99
  #   TMLE[["TMLE"]][[1]][["SE.TMLE"]]
  #   # [1] 0.009949874 0.009949874 0.009949874

  #   TMLE[["TMLE"]][[2]][["St.TMLE"]]
  #   # [1] 0.9900000 0.9789573 0.9679731
  #   TMLE[["TMLE"]][[2]][["SE.TMLE"]]
  #   # [1] 0.009949874 0.015098587 0.015938640
  # })
