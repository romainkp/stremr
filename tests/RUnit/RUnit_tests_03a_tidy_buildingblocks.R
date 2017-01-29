test.GRID.h2o.xgboost.10Kdata <- function() {
  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  if (reqxgb) {
    `%+%` <- function(a, b) paste0(a, b)
    # options(stremr.verbose = TRUE)
    # options(gridisl.verbose = FALSE)
    options(stremr.verbose = FALSE)
    options(gridisl.verbose = FALSE)
    library("data.table")
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

    params_g <- gridisl::defModel(estimator = "xgboost__glm",
                                  family = "binomial",
                                  nrounds = 10,
                                  early_stopping_rounds = 2)
                  # gridisl::defModel(estimator = "h2o__glm",
                  #                   family = "binomial",
                  #                   lambda_search = TRUE,
                  #                   param_grid = list(
                  #                     alpha = c(0.5)
                  #                   )) +
                  # gridisl::defModel(estimator = "xgboost__gbm",
                  #               family = "binomial",
                  #               search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
                  #               seed = 23,
                  #               # learning_rate = 0.1, # learning_rate = 0.01, # learning_rate = 0.05,
                  #               # nthread = 1,
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
                  #                   ## Maximum delta step we allow each tree’s weight estimation to be.
                  #                   ## If the value is set to 0, it means there is no constraint.
                  #                   ## If it is set to a positive value, it can help making the update step more conservative.
                  #                   ## Might help in logistic regression when class is extremely imbalanced.
                  #                   ## Set it to value of 1-10 to help control the update
                  #                   max_delta_step = c(0, 1, 2, 5, 10)
                  #                   )
                  #               )

    OData <- fitPropensity(OData,
                            gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            models_CENS = params_g, models_TRT = params_g, models_MONITOR = params_g,
                            fit_method = "cv",
                            fold_column = "fold_ID")

  ## ----------------------------------------------------------------------------------------------------
  ## Analysis by intervention
  ## **** makes it easier to read the individual analysis ****
  ## ----------------------------------------------------------------------------------------------------
  require("magrittr")
  # require("tidyverse")
  # require("tidyr")
  # require("dplyr")
  # require("purrr")
  library("ggplot2")
  library("tibble")
  library("tidyr")
  library("readr")
  library("purrr")
  library("dtplyr")
  require("dplyr")

  # t_periods <- 0:8
  t_periods <- 0:8
  t_breaks = c(1:8,12,16)-1
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t_periods)+1))

  ## ------------------------------------------------------------
  ## **** 1) As a first step define a grid of all possible parameter combinations (for all estimators)
  ## **** 2) This dataset is to be saved and will be later merged in with all analysis
  ## ------------------------------------------------------------
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh", "gTI.dlow"),
                  trunc_wt = c(FALSE, TRUE),
                  stratifyQ_by_rule = c(TRUE, FALSE)) %>%
                  cross_d() %>%
                  arrange(stratifyQ_by_rule) %>%
                  mutate(trunc_MSM = map_dbl(trunc_wt, ~ ifelse(.x, 50, Inf))) %>%
                  mutate(trunc_TMLE = trunc_MSM*10)

  ## ------------------------------------------------------------
  ## IPW ANALYSIS
  ## **** 3) For each individual analysis do filter()/subset()/etc to create a grid of parameters specific to given estimator
  ## ------------------------------------------------------------
  IPW <-  analysis %>%
          rename(trunc_weight = trunc_MSM) %>%
          distinct(intervened_TRT, trunc_weight) %>%

          group_by(intervened_TRT) %>%
          mutate(wts_data = map(first(intervened_TRT), getIPWeights, OData = OData)) %>%

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
                          t_breaks = t_breaks,
                          use_weights = FALSE,
                          glm_package = "speedglm"))) %>%
          mutate(MSM.crude = map(MSM.crude, "estimates")) %>%

          ## IPW-MSM for hazard (smoothing over time-intervals in t_breaks):
          mutate(MSM = map2(wts_data, trunc_weight,
              ~ survMSM(wts_data = .x,
                          trunc_weights = .y,
                          OData = OData,
                          t_breaks = t_breaks,
                          glm_package = "speedglm"))) %>%
          mutate(MSM = map(MSM, "estimates")) %>%
          rename(trunc_MSM = trunc_weight)

  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  GCOMP <-analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule) %>%

          mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
              ~ fitSeqGcomp(intervened_TRT = .x,
                            stratifyQ_by_rule = .y,
                            t_periods = t_periods,
                            OData = OData,
                            Qforms = Qforms)
                          )) %>%
          mutate(GCOMP = map(GCOMP, "estimates"))

  # GCOMP <-GCOMP %>%
  #         mutate(plotGCOMP = map(GCOMP, ~ ggsurv(estimates = .x))) %>%

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  TMLE <- analysis %>%
          rename(trunc_weight = trunc_TMLE) %>%
          distinct(intervened_TRT, stratifyQ_by_rule, trunc_weight)

  TMLE <- TMLE %>%
          mutate(TMLE = pmap(TMLE, fitTMLE,
                             t_periods = t_periods,
                             OData = OData,
                             Qforms = Qforms)) %>%
          mutate(TMLE = map(TMLE, "estimates")) %>%
          rename(trunc_TMLE = trunc_weight)

  ## ------------------------------------------------------------
  ## COMBINE ALL ANALYSES INTO A SINGLE DATASET
  ## ------------------------------------------------------------
  results <-  analysis %>%
              left_join(IPW) %>%
              left_join(GCOMP) %>%
              left_join(TMLE)

  ## ------------------------------------------------------------
  ## REMOVE wts_data (no longer needed)
  ## TO DO: NEED TO CREATE SUMMARY TABLES FIRST AND REPLACE wts_data
  ##        get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE)
  ## Nest each estimator by treatment regimen (we now only show the main analysis rows)
  ## ------------------------------------------------------------
  results_nowts <- results %>%
              select(-wts_data) %>%
              nest(intervened_TRT, NPMSM, MSM.crude, MSM, GCOMP, TMLE, .key = "estimates")

  # results_nowts[["estimates"]]
  # results_nowts[1, ][["estimates"]][[1]][["TMLE"]]

  ## ------------------------------------------------------------
  ## Calculate RDs (contrasting all interventions, for each analysis row & estimator)
  ## The RDs data no longer needs the intervened_TRT column
  ## ------------------------------------------------------------
  fin_results <-  results_nowts %>%
                  mutate(RDs =
                      map(estimates,
                          ~ select(.x, -intervened_TRT) %>%
                              map(~ get_RDs(.x)) %>%
                              as_tibble()
                          )
                      )

  ## clean up by removing the subject-level IC estimates for EVERY SINGLE ESTIMATE / ANALYSIS
  ## THIS IS A SIDE-EFFECT FUNCTION
  res <- fin_results[["estimates"]] %>%
      map(
          ~ select(.x, -intervened_TRT) %>%
              map(
                  ~ map(.x, ~ suppressWarnings(.x[, ("IC.St") := NULL])))
          )

      ## equivalent to above but with explicitely defined inside function::
      # fin_results_2 <- results_nowts %>%
      #     mutate(RDs = map(estimates, function(.df) {
      #         res <- select(.df, -intervened_TRT) %>%
      #         map(~ get_RDs(.x))
      #         browser()
      #         as_tibble(res)
      #         return(res)
      #         })
      #     )
      # fin_results[["RDs"]][[1]][["MSM"]]
      # fin_results_2[["estimates"]][[1]]
      # fin_results_2[["RDs"]]

      ## equivalent to above two
      # fin_results <-results_nowts %>%
      #             as_tibble() %>%
      #             group_by(trunc_wt, stratifyQ_by_rule) %>%
      #             mutate(NPMSM = get_RDs(NPMSM, "St.NPMSM", getSEs = FALSE)) %>%
      #             mutate(MSM.crude = get_RDs(MSM.crude, "St.MSM", getSEs = TRUE)) %>%
      #             mutate(MSM = get_RDs(MSM, "St.MSM", getSEs = TRUE)) %>%
      #             mutate(TMLE = get_RDs(TMLE, "St.TMLE", getSEs = FALSE)) %>%
      #             ungroup() %>%
      #             distinct(trunc_wt, stratifyQ_by_rule, NPMSM, MSM.crude, MSM, TMLE)

      # results_nowts[["estimates"]][[1]][["TMLE"]]
      # fin_results[["estimates"]][[1]][["TMLE"]]
      ## equivalent to above:
      # fin_results[["estimates"]] %>%
      #     map( function(.df) {
      #         res <- select(.df, -intervened_TRT) %>%
      #         map( ~ map(.x, ~ .x[, ("IC.St") := NULL]))
      #         browser()
      #         return(res)
      #         })

  ## ------------------------------------------------------------
  ## PLOTTING RDs
  ## ------------------------------------------------------------
      estimator <- "TMLE"
      RDplot <-   fin_results[["RDs"]][[1]][[estimator]][[1]] %>%
                  ggRD %>%
                  print

      ## across scenarios / estimators
      longRDs <- fin_results %>%
          select(trunc_wt, stratifyQ_by_rule, trunc_MSM, trunc_TMLE, RDs) %>%
          unnest(RDs) %>%
          gather(key = est, value = RDs, NPMSM, MSM.crude, MSM, GCOMP, TMLE) %>%
          filter(est %in% c("TMLE", "MSM")) %>%
          mutate(RDplot = map(RDs, ~ ggRD(.x)))

  # longRDs[["RDplot"]]

  longRDs %>%
      trelliscopejs::trelliscope(name = "test", panel_col = "RDplot")
      longRDs







  ## NO LONGER USED
  # ## Calculate RDs and remove the subject-level IC estimates:
  # IPW <- IPW %>%
  #         group_by(trunc_MSM) %>%
  #         mutate(MSM.RD = get_RDs(MSM, "St.MSM", getSEs = TRUE)) %>%
  #         ungroup %>%
  #         mutate(MSM = map(MSM, ~ .x[, ("IC.St") := NULL]))
  # ## Calculate RDs and remove the subject-level IC estimates:
  # GCOMP <- GCOMP %>%
  #         group_by(stratifyQ_by_rule) %>%
  #         mutate(GCOMP.RD = get_RDs(GCOMP, "St.GCOMP", getSEs = FALSE)) %>%
  #         ungroup
  # ## Calculate RDs and remove the subject-level IC estimates:
  # TMLE <- TMLE %>%
  #         group_by(stratifyQ_by_rule, trunc_TMLE) %>%
  #         mutate(TMLE.RD = get_RDs(TMLE, "St.TMLE", getSEs = TRUE)) %>%
  #         ungroup %>%
  #         mutate(TMLE = map(TMLE, ~ .x[, ("IC.St") := NULL]))
  # resultsRDs <- results %>%
  #     group_by(trunc_wt, stratifyQ_by_rule) %>%
  #     mutate(TMLE.RD = get_RDs(TMLE, "St.TMLE", getSEs = TRUE)) %>%
  #     ungroup



  results[["TMLE"]]
  str(lapply(results[["TMLE"]], attributes))
  ggsurv(results[["TMLE"]])

  results[["t_breaks"]]
  results[["NPMSM"]]
  results[["MSM"]]
  results[["GCOMP"]]
  results[["TMLE"]]
  setDT(results)
  attributes(results[["MSM"]][[1]])

  MSM.RDtables2 = get_RDs(results[["MSM"]], "St.MSM", getSEs = TRUE)
  TMLE.RDtables2 = get_RDs(results[["TMLE"]], "St.TMLE", getSEs = TRUE)

      pl <- ggsurv(estimates = results[["NPMSM"]])
      pl2 <- ggsurv(estimates = results[["MSM"]])
      pl3 <- ggsurv(estimates = results[["GCOMP"]])
      pl4 <- ggsurv(estimates = results[["TMLE"]])
      pl4 <- ggsurv(estimates = results[["TMLE"]][[1]])
}
}