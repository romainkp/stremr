context("origami Super Learner")

  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  reqh2o <- requireNamespace("h2o", quietly = TRUE)
  if (!reqxgb || !reqh2o) return(TRUE)

  ## -----------------------------------------------------------------------
  ## Analyses by intervention
  ## **** makes it easier to read the individual analyses ****
  ## -----------------------------------------------------------------------
  `%+%` <- function(a, b) paste0(a, b)
  library("stremr")
  library("data.table")
  library("magrittr")
  library("ggplot2")
  library("tibble")
  library("tidyr")
  library("purrr")
  library("dplyr")
  library("sl3")

  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  options(sl3.verbose = FALSE)
  options(condensier.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = TRUE)
  # options(sl3.verbose = TRUE)
  # options(condensier.verbose = TRUE)
  
  data.table::setDTthreads(1)

  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # select only the first 100 IDs
  Odat_DT <- Odat_DT[ID %in% (1:100), ]
  setkeyv(Odat_DT, cols = c("ID", "t"))

  ## -----------------------------------------------------------------------
  ## Define some summaries (lags C[t-1], A[t-1], N[t-1])
  ## -----------------------------------------------------------------------
  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- lagnodes %+% ".tminus1"
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

  gform_CENS <- c("C ~ highA1c + t")
  gform_MONITOR <- "N ~ 1"

  ## ------------------------------------------------------------
  ## **** As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** This dataset is to be saved and will be later merged in with all analysis
  ## ------------------------------------------------------------
  trunc_IPW <- 10
  # tvals <- 0:6
  tvals <- 0:2
  tmax <- 13
  ## number of folds for CV:
  # nfolds <- 6
  nfolds <- 4
  tbreaks = c(1:8,11,14)-1

  ## This dataset defines all parameters that we like to vary in this analysis (including different interventions)
  ## That is, each row of this dataset corresponds with a single analysis, for one intervention of interest.
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh"),
                  trunc_wt = c(FALSE, TRUE),
                  stratifyQ_by_rule = c(TRUE, FALSE)) %>%
                  cross_df() %>%
                  arrange(stratifyQ_by_rule) %>%
                  mutate(nfolds = as.integer(nfolds)) %>%
                  mutate(trunc_MSM = map_dbl(trunc_wt, ~ ifelse(.x, trunc_IPW, Inf))) %>%
                  mutate(trunc_TMLE = trunc_MSM*10)

  ## ----------------------------------------------------------------
  ## IMPORT DATA
  ## ----------------------------------------------------------------
  # library("h2o")
  # h2o::h2o.init(nthreads = 2)

  OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  OData <- define_CVfolds(OData, nfolds = nfolds, fold_column = "fold_ID", seed = 12345)
  OData$dat.sVar[]
  OData$fold_column <- NULL
  OData$nfolds <- NULL
  fold_column <- "fold_ID"

  ## ----------------------------------------------------------------
  ## Define ensemble of models for fitting propensity scores (g).
  ## ----------------------------------------------------------------
  ## This example uses a discrete SuperLearner: best model will be selected on the basis of CV-MSE.
  ## Use cross-validation to select best model for g (set 'fit_method_g <- "none"' to just fit a single model w/out CV)
  # fit_method_g <- "origamiSL"
  fit_method_g <- "cv"
  # fit_method_g <- "none"

  ## Note that 'interactions' CANNOT be used with h2o (for now).
  ## The only learners that allow interactions are: "glm" ,"speedglm", "xgboost".
  models_g <-
     # defModel(estimator = "h2o__glm", family = "binomial",
     #                   # lambda_search = FALSE,
     #                   # nlambdas = 5,
     #                   param_grid = list(
     #                    alpha = 0
     #                    # alpha = c(0.5)
     #                   ))
    # defModel(estimator = "xgboost__glm", family = "binomial", nthread = 1) +
    # defModel(estimator = "speedglm__glm", family = "quasibinomial") +
    # defModel(estimator = "xgboost__gbm", family = "quasibinomial", nrounds = 20)
    Lrnr_sl$new(learners = Stack$new(
                Lrnr_xgboost$new(nthread = 1, nrounds = 20),
                Lrnr_glm_fast$new()
                ),
                metalearner = Lrnr_solnp$new()
            )

  ## ------------------------------------------------------------------------
  ## Define models for iterative G-COMP (Q) -- PARAMETRIC LOGISTIC REGRESSION
  ## ------------------------------------------------------------------------
  ## regression formulas for Q's:
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(tvals)+1))

  ## no cross-validation model selection, just fit a single model specified below
  # fit_method_Q <- "origamiSL"
  fit_method_Q <- "cv"

  ## Use speedglm to fit all Q.
  ## NOTE: it is currently not possible to use fit_method_Q <- "cv" with speedglm or glm.
  ## To perform cross-validation with GLM use 'estimator="h2o__glm"' or 'estimator="xgboost__glm"'
  models_Q <-
    Lrnr_sl$new(learners = Stack$new(
                Lrnr_xgboost$new(nthread = 1, nrounds = 20),
                Lrnr_glm_fast$new()
                ),
                metalearner = Lrnr_solnp$new()
            )

  ## ------------------------------------------------------------------------
  ## Alternative specifications of Q models for iterative G-COMP (Q) -- NONPARAMETRIC REGRESSION (GBM) + REGULARIZED GLM
  ## ------------------------------------------------------------------------
  # ## This example uses a discrete SuperLearner: best model will be selected on the basis of CV-MSE.
  # fit_method_Q <- "cv"
  # models_Q <<-
  #   defModel(estimator = "xgboost__gbm", family = "binomial",
  #                       nrounds = 200,
  #                       early_stopping_rounds = 3,
  #                       interactions = list(c("TI", "highA1c"), c("TI", "CVD"), c("TI", "lastNat1")),
  #                       learning_rate = .1,
  #                       max_depth = 5,
  #                       gamma = .5,
  #                       colsample_bytree = 0.8,
  #                       subsample = 0.8,
  #                       lambda = 2,
  #                       alpha = 0.5,
  #                       max_delta_step = 2) +
  #   defModel(estimator = "xgboost__glm",
  #                     family = "quasibinomial",
  #                     seed = 23,
  #                     nrounds = 500,
  #                     early_stopping_rounds = 5,
  #                     interactions = list(c("TI", "highA1c"), c("TI", "CVD"), c("TI", "lastNat1")),
  #                     param_grid = list(
  #                       alpha = c(.0, 0.5, 1.0),
  #                       lambda = c(.01, .1, .5, .9, 1.5, 5)
  #                       )
  #                     )

  ## ----------------------------------------------------------------
  ## Fit propensity score models.
  ## We are using the same model ensemble defined in models_g for censoring, treatment and monitoring mechanisms.
  ## ----------------------------------------------------------------
  # test_that("fitting g w/ Super Learner", {
  OData <- fitPropensity(OData,
                          gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT,
                          models_CENS = models_g, models_TRT = models_g,
                          fit_method = fit_method_g,
                          fold_column = fold_column)
  # })

  ## ------------------------------------------------------------
  ## RUN IPW ANALYSES
  ## **** For each individual analysis do filter()/subset()/etc to create a grid of parameters specific to given estimator
  ## ------------------------------------------------------------
# test_that("evaluating IPW", {
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
# })

  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  # test_that("fitting Q w/ origami Super Learner, with byfold = TRUE", {
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
                        byfold = FALSE))) %>%
        mutate(GCOMP = map(GCOMP, "estimates"))
  # })

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  # TMLE <- CVTMLE <- analysis %>%
  #         rename(trunc_weight = trunc_TMLE) %>%
  #         distinct(intervened_TRT, stratifyQ_by_rule, trunc_weight)

  # TMLE_time <- system.time({
  #   TMLE <- TMLE %>%
  #         mutate(TMLE = pmap(TMLE, fit_TMLE,
  #                            tvals = tvals,
  #                            OData = OData,
  #                            models = models_Q,
  #                            Qforms = Qforms,
  #                            fit_method = fit_method_Q,
  #                            fold_column = fold_column)) %>%
  #         mutate(TMLE = map(TMLE, "estimates")) %>%
  #         rename(trunc_TMLE = trunc_weight)
  # })
  # TMLE_time_hrs <- TMLE_time[3]/60/60

  # CVTMLE_time <- system.time({
  #   CVTMLE <- CVTMLE %>%
  #         mutate(CVTMLE = pmap(CVTMLE, fit_CVTMLE,
  #                            tvals = tvals,
  #                            OData = OData,
  #                            models = models_Q,
  #                            Qforms = Qforms,
  #                            fit_method = fit_method_Q,
  #                            fold_column = fold_column)) %>%
  #         mutate(CVTMLE = map(CVTMLE, "estimates")) %>%
  #         rename(trunc_TMLE = trunc_weight)
  # })
  # CVTMLE_time_hrs <- TMLE_time[3]/60/60


  ## ------------------------------------------------------------
  ## COMBINE ALL ANALYSES INTO A SINGLE DATASET
  ## ------------------------------------------------------------
  results <-  analysis %>%
              left_join(IPW) %>%
              left_join(GCOMP)
              # left_join(TMLE) %>%
              # left_join(CVTMLE)

  ## Nest each estimator by treatment regimen (we now only show the main analysis rows)
  results <- results %>%
              nest(intervened_TRT, NPMSM, MSM.crude, MSM, GCOMP, .key = "estimates")
              # nest(intervened_TRT, NPMSM, MSM.crude, MSM, GCOMP, TMLE, CVTMLE, .key = "estimates")

  ## Calculate RDs (contrasting all interventions, for each analysis row & estimator).
  ## The RDs data no longer needs the intervened_TRT column
  results <-  results %>%
              mutate(RDs =
                map(estimates,
                  ~ select(.x, -intervened_TRT) %>%
                  map(~ get_RDs(.x)) %>%
                  as_tibble()
                  ))

  ## ------------------------------------------------------------
  ## Uncomment and run this to remove the individual EIC estimates from MSM and TMLE estimates.
  ## This is useful if saving results to a file. The EIC take up a lot of space, hence removing them
  ## will significantly reduce the final file size.
  ## ------------------------------------------------------------
  ## Clean up by removing the subject-level IC estimates for EVERY SINGLE ESTIMATE / ANALYSIS
  ## WARNING: THIS IS A SIDE-EFFECT FUNCTION!
  # res <- results[["estimates"]] %>%
  #     map(
  #       ~ select(.x, -intervened_TRT) %>%
  #       map(
  #         ~ map(.x,
  #           ~ suppressWarnings(.x[, ("IC.St") := NULL]))))
  # rm(res)

  ## ------------------------------------------------------------
  ## Add models used for g and Q. Create the final analysis file.
  ## Add IPWtabs
  ## ------------------------------------------------------------

  # results <- results %>%
  #            left_join(IPWtabs) %>%
  #            mutate(fit_method_g = fit_method_g) %>%
  #            mutate(fit_method_Q = fit_method_Q) %>%
  #            mutate(models_g = map(fit_method_g, ~ models_g)) %>%
  #            mutate(models_Q = map(fit_method_Q, ~ models_Q)) %>%
  #            mutate(run_time = map(trunc_wt,
  #             ~ tibble(IPW_time_hrs = IPW_time_hrs, GCOMP_time_hrs = GCOMP_time_hrs, TMLE_time_hrs = TMLE_time_hrs)))

  ## ------------------------------------------------------------
  ## VARIOUS WAYS OF PLOTTING SURVIVAL CURVES
  ## ------------------------------------------------------------
  ests <- "GCOMP"
  SURVplot <- results[1, ][["estimates"]][[1]][[ests]] %>%
              ggsurv %>%
              print

  # results[["estimates"]]

  # GCOMP <-GCOMP %>%
  #         mutate(plotGCOMP = map(GCOMP, ~ ggsurv(estimates = .x))) %>%

  ## FLATTEN THE results data (long format)
  ## then add a new column with ggplot survival plot objects for THESE 3 estimators;
  # ests <- c("MSM", "GCOMP", "TMLE")
  # longSURV <- results %>%
  #             select(trunc_wt, stratifyQ_by_rule, trunc_MSM, trunc_TMLE, estimates) %>%
  #             unnest(estimates) %>%
  #             gather(key = est, value = estimates, NPMSM, MSM.crude, MSM, GCOMP, TMLE, CVTMLE) %>%
  #             filter(est %in% ests) %>%
  #             select(-intervened_TRT) %>%
  #             nest(estimates, .key = "estimates") %>%
  #             mutate(SURVplot = map(estimates, ~ ggsurv(.x)))

  ## Visualize all survival curves at once with a single interactive trelliscope panel
  ## (install via: devtools::install_github("hafen/trelliscopejs"))
  # reqtrell <- requireNamespace("trelliscopejs", quietly = TRUE)
  # if (reqtrell && FALSE) {
  #   longSURV %>%
  #     trelliscopejs::trelliscope(name = "test", panel_col = "SURVplot")
  #     longSURV
  # }

  ## ------------------------------------------------------------
  ## VARIOUS WAYS OF PLOTTING RD TABLEs
  ## ------------------------------------------------------------
  ## THE TABLE OF RDs FOR TMLE:
  # results %>% filter(trunc_wt == TRUE, stratifyQ_by_rule == TRUE) %>% select(RDs) %>% unnest(RDs) %>% select(TMLE) %>% unnest(TMLE)

  # ## PLOT RD FOR TMLE FOR SPECIFIC SCENARIO:
  # ests <- "TMLE"
  # RDplot <-   results[["RDs"]][[1]][[ests]][[1]] %>%
  #             ggRD(t_int_sel = 1:5) %>%
  #             print

  # ## GENERERATE RD PLOTS across all scenarios for these two estimators:
  # ests <- c("MSM", "TMLE")
  # longRDs <- results %>%
  #             select(trunc_wt, stratifyQ_by_rule, trunc_MSM, trunc_TMLE, RDs) %>%
  #             unnest(RDs) %>%
  #             gather(key = est, value = RDs, NPMSM, MSM.crude, MSM, GCOMP, TMLE, CVTMLE) %>%
  #             filter(est %in% ests) %>%
  #             mutate(RDplot = map(RDs, ~ ggRD(.x)))

  # reqtrell <- requireNamespace("trelliscopejs", quietly = TRUE)
  # if (reqtrell && FALSE) {
  #   longRDs %>%
  #     trelliscopejs::trelliscope(name = "test", panel_col = "RDplot")
  #     longRDs
  # }

  # h2o::h2o.shutdown(prompt = FALSE)
# }