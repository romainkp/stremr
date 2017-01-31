test.GRID.h2o.xgboost.10Kdata <- function() {
  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  if (!reqxgb) return()

  ## ----------------------------------------------------------------------------------------------------
  ## Analysis by intervention
  ## **** makes it easier to read the individual analysis ****
  ## ----------------------------------------------------------------------------------------------------
  `%+%` <- function(a, b) paste0(a, b)
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = FALSE)
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)

  library("data.table")
  library("magrittr")
  library("ggplot2")
  library("tibble")
  library("tidyr")
  library("purrr")
  library("dtplyr")
  library("dplyr")

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

  fold_column <- "fold_ID"
  fit_method_g <- "cv"
  models_g <- gridisl::defModel(estimator = "xgboost__glm",
                                family = "binomial",
                                nrounds = 10,
                                early_stopping_rounds = 2)
                # gridisl::defModel(estimator = "xgboost__gbm",
                #               family = "binomial",
                #               search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
                #               seed = 23,
                #               nrounds = 200, early_stopping_rounds = 10,
                #               param_grid = list(
                #                   learning_rate = c(.1, .3, .5), # .05,
                #                   max_depth = c(seq(3, 19, 4), 25),
                #                   min_child_weight = c(1, 3, 5, 7),
                #                   gamma = c(.0, .05, seq(.1, .9, by=.2), 1),
                #                   colsample_bytree = c(.6, .8, 1), # .4,
                #                   subsample = c(.5, .75, 1),
                #                   lambda = c(.1, .5, 2, 5), # lambda = c(1,2,5),
                #                   alpha = c(0, .1, .5),
                #                   max_delta_step = c(0, 1, 2, 5, 10)
                #                   )
                #               )

  fit_method_Q <- "none"
  models_Q <-  gridisl::defModel(estimator = "speedglm__glm", family = "quasibinomial")

  OData <- fitPropensity(OData,
                          gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                          models_CENS = models_g, models_TRT = models_g, models_MONITOR = models_g,
                          fit_method = fit_method_g,
                          fold_column = fold_column)

  trunc_IPW <- 10
  tvals <- 0:8
  tmax <- 13
  nfolds <- 10 ## number of folds for CV
  # tbreaks = c(1:8,12,16)-1
  tbreaks = c(1:8,11,14)-1
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(tvals)+1))

  ## ------------------------------------------------------------
  ## **** As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** This dataset is to be saved and will be later merged in with all analysis
  ## ------------------------------------------------------------
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh", "gTI.dlow"),
                  trunc_wt = c(FALSE, TRUE),
                  stratifyQ_by_rule = c(TRUE, FALSE)) %>%
                  cross_d() %>%
                  arrange(stratifyQ_by_rule) %>%
                  mutate(nfolds = as.integer(nfolds)) %>%
                  mutate(trunc_MSM = map_dbl(trunc_wt, ~ ifelse(.x, trunc_IPW, Inf))) %>%
                  mutate(trunc_TMLE = trunc_MSM*10)

  ## ------------------------------------------------------------
  ## IPW ANALYSIS
  ## **** For each individual analysis do filter()/subset()/etc to create a grid of parameters specific to given estimator
  ## ------------------------------------------------------------
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
          rename(trunc_MSM = trunc_weight) %>%
          select(-wts_data)

  ## save IPW tables (will be later merged with main results dataset)
  IPWtabs <-  analysis %>%
               left_join(IPW) %>%
               distinct(intervened_TRT, trunc_MSM, wts_tabs, FUPtimes_tabs) %>%
               nest(intervened_TRT, wts_tabs, FUPtimes_tabs, .key = "IPWtabs")

  IPW <- IPW %>% select(-wts_tabs, -FUPtimes_tabs)

  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  GCOMP <-analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule) %>%
          mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
            ~ fitSeqGcomp(intervened_TRT = .x,
                          stratifyQ_by_rule = .y,
                          tvals = tvals,
                          OData = OData,
                          models = models_Q,
                          Qforms = Qforms,
                          fit_method = fit_method_Q,
                          fold_column = fold_column))) %>%
          mutate(GCOMP = map(GCOMP, "estimates"))

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  TMLE <- analysis %>%
          rename(trunc_weight = trunc_TMLE) %>%
          distinct(intervened_TRT, stratifyQ_by_rule, trunc_weight)

  TMLE <- TMLE %>%
          mutate(TMLE = pmap(TMLE, fitTMLE,
                             tvals = tvals,
                             OData = OData,
                             models = models_Q,
                             Qforms = Qforms,
                             fit_method = fit_method_Q,
                             fold_column = fold_column)) %>%
          mutate(TMLE = map(TMLE, "estimates")) %>%
          rename(trunc_TMLE = trunc_weight)

  ## ------------------------------------------------------------
  ## COMBINE ALL ANALYSES INTO A SINGLE DATASET
  ## ------------------------------------------------------------
  results <-  analysis %>%
              left_join(IPW) %>%
              left_join(GCOMP) %>%
              left_join(TMLE)

  ## Nest each estimator by treatment regimen (we now only show the main analysis rows)
  results <- results %>%
              # select(-wts_tabs, -FUPtimes_tabs) %>%
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
  results <- results %>%
             mutate(fit_method_g = fit_method_g) %>%
             mutate(fit_method_Q = fit_method_Q) %>%
             mutate(models_g = map(fit_method_g, ~ I(models_g))) %>%
             mutate(models_Q = map(fit_method_Q, ~ I(models_Q))) %>%
             left_join(IPWtabs)

  ## ------------------------------------------------------------
  ## PLOTTING SURVIVAL CURVES
  ## ------------------------------------------------------------
  ests <- "TMLE"
  SURVplot <- results[1, ][["estimates"]][[1]][[ests]] %>%
              ggsurv %>%
              print

  results[["estimates"]]

  # GCOMP <-GCOMP %>%
  #         mutate(plotGCOMP = map(GCOMP, ~ ggsurv(estimates = .x))) %>%
  ests <- c("MSM", "TMLE")
  longSURV <- results %>%
              select(trunc_wt, stratifyQ_by_rule, trunc_MSM, trunc_TMLE, estimates) %>%
              unnest(estimates) %>%
              gather(key = est, value = estimates, NPMSM, MSM.crude, MSM, GCOMP, TMLE) %>%
              filter(est %in% ests) %>%
              select(-intervened_TRT) %>%
              nest(estimates, .key = "estimates") %>%
              mutate(SURVplot = map(estimates, ~ ggsurv(.x)))

  reqtrell <- requireNamespace("trelliscopejs", quietly = TRUE)
  if (reqtrell) {
    longSURV %>%
      trelliscopejs::trelliscope(name = "test", panel_col = "SURVplot")
      longSURV
  }

  ## ------------------------------------------------------------
  ## PLOTTING RDs
  ## ------------------------------------------------------------
  ests <- "TMLE"
  RDplot <-   results[["RDs"]][[1]][[ests]][[1]] %>%
              ggRD(t_int_sel = 1:5) %>%
              print

  ## across all scenarios for two estimators
  ests <- c("MSM", "TMLE")
  longRDs <- results %>%
              select(trunc_wt, stratifyQ_by_rule, trunc_MSM, trunc_TMLE, RDs) %>%
              unnest(RDs) %>%
              gather(key = est, value = RDs, NPMSM, MSM.crude, MSM, GCOMP, TMLE) %>%
              filter(est %in% ests) %>%
              mutate(RDplot = map(RDs, ~ ggRD(.x)))

  reqtrell <- requireNamespace("trelliscopejs", quietly = TRUE)
  if (reqtrell) {
    longRDs %>%
      trelliscopejs::trelliscope(name = "test", panel_col = "RDplot")
      longRDs
  }

}