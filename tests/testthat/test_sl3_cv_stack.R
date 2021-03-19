context("Testing sl3 cross-validation and super learner")
  # devtools::install_github("jeremyrcoyle/sl3")
  library("stremr")
  library("sl3")
  library("SuperLearner")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  options(sl3.verbose = FALSE)
  library("data.table")
  library("magrittr")
  library("ggplot2")
  # library("tibble")
  # library("tidyr")
  library("purrr")
  library("dplyr")

  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # select only the first 100 IDs
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

  ## remove indicators of censoring and monitoring events:
  Odat_DT[, "N" := NULL]
  # Odat_DT[, "C" := NULL]

  ## ------------------------------------------------------------------
  ## Propensity score models for Treatment, Censoring & Monitoring
  ## ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                      # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))
  ## ------------------------------------------------------------------------
  ## Outcome models for iterative G-COMP (Q)
  ## ------------------------------------------------------------------------
  tvals <- 2
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + lastNat1 + TI + TI.tminus1", (max(tvals)+1))

  ## ------------------------------------------------------------
  ## **** As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** This dataset is to be saved and will be later merged in with all analysis
  ## ------------------------------------------------------------
  ## This dataset defines all parameters that we like to vary in this analysis (including different interventions)
  ## That is, each row of this dataset corresponds with a single analysis, for one intervention of interest.
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh"), stratifyQ_by_rule = c(FALSE)) %>%
                  cross_df() %>%
                  arrange(stratifyQ_by_rule)

  ## ----------------------s--------------------------------------------------
  ## Define models for fitting propensity scores (g) and outcome model (Q)
  ## ------------------------------------------------------------------------
  lrn_glm_base <- Lrnr_glm$new(family = binomial())
  lrn_glm <- Lrnr_glm_fast$new()
  lrn_glm_sm <- Lrnr_glm_fast$new(covariates = c("CVD"))
  # lrn_glmnet_binom <- Lrnr_pkg_SuperLearner$new("SL.glmnet", family = binomial())
  # lrn_glmnet_binom <- Lrnr_glmnet$new(family = "binomial", nlambda = 5)
  # lrn_glmnet_gaus <- Lrnr_pkg_SuperLearner$new("SL.glmnet", family = "gaussian")
  # lrn_glmnet_gaus <- Lrnr_glmnet$new(family = "gaussian", nlambda = 5)
  sl_g <- Lrnr_sl$new(learners = Stack$new(lrn_glm_base, lrn_glm, lrn_glm_sm),
                    metalearner = Lrnr_nnls$new())
  sl_Q <- Lrnr_sl$new(learners = Stack$new(lrn_glm, lrn_glm_sm),
                    metalearner = Lrnr_nnls$new())

  fit_method_g <- "cv"
  # fit_method_g <- "none"    
  ## no cross-validation model selection, just fit a single model specified below
  fit_method_Q <- "none"

  ## ----------------------------------------------------------------
  ## Fit propensity score models.
  ## We are using the same model ensemble defined in sl_g for censoring, treatment and monitoring mechanisms.
  ## ----------------------------------------------------------------
  OData <- stremr::importData(Odat_DT, ID = "ID", t_name = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), TRT = "TI", OUTCOME = "Y.tplus1") %>% 
           stremr::define_CVfolds(nfolds = 4, fold_column = "fold_ID")

  OData <- fitPropensity(OData,
                         gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT,
                         models_TRT = sl_g,
                         fit_method = fit_method_g
                        )

  ## Get the dataset with weights:
  wts_data <- getIPWeights(intervened_TRT = "gTI.dlow", OData = OData)

  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  # GCOMP <-analysis %>%
  #       distinct(intervened_TRT, stratifyQ_by_rule) %>%
  #       mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
  #         ~ fit_GCOMP(intervened_TRT = .x,
  #                       stratifyQ_by_rule = .y,
  #                       tvals = tvals,
  #                       OData = OData,
  #                       models = models_Q,
  #                       Qforms = Qforms,
  #                       fit_method = fit_method_Q
  #                       ))) %>%
  #       mutate(GCOMP = map(GCOMP, "estimates"))


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
  TMLE <- analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule)

  TMLE <- TMLE %>%
        mutate(TMLE = pmap(TMLE, fit_TMLE,
                           tvals = tvals,
                           OData = OData,
                           models = sl_Q,
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
