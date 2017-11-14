context("Fitting with no Monitoring and / or no Censoring indicators")

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

  ## remove indicators of censoring and monitoring events:
  Odat_DT[, "N" := NULL]
  Odat_DT[, "C" := NULL]

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
  tvals <- 0:2
  tmax <- 13
  tbreaks = c(1:8,11,14)-1

  ## This dataset defines all parameters that we like to vary in this analysis (including different interventions)
  ## That is, each row of this dataset corresponds with a single analysis, for one intervention of interest.
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh"),
                  stratifyQ_by_rule = c(TRUE)) %>%
                  cross_df() %>%
                  arrange(stratifyQ_by_rule)

  ## ------------------------------------------------------------------------
  ## Define models for fitting propensity scores (g) -- PARAMETRIC LOGISTIC REGRESSION
  ## ------------------------------------------------------------------------
  fit_method_g <- "none"
  # models_g <- defModel(estimator = "speedglm__glm", family = "quasibinomial")

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

  ## ----------------------------------------------------------------
  ## Fit propensity score models.
  ## We are using the same model ensemble defined in models_g for censoring, treatment and monitoring mechanisms.
  ## ----------------------------------------------------------------

  OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), TRT = "TI", OUTCOME = outcome)
  OData <- fitPropensity(OData,
                          gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT,
                          # models_TRT = models_g,
                          fit_method = fit_method_g
                          )

  ## Get the dataset with weights:
  wts_data <- getIPWeights(intervened_TRT = "gTI.dlow", OData = OData, tmax = tmax)

  ## ------------------------------------------------------------
  ## RUN IPW ANALYSES
  ## **** For each individual analysis do filter()/subset()/etc to create a grid of parameters specific to given estimator
  ## ------------------------------------------------------------
  # test_that("evaluating IPW", {
    IPW <-  analysis %>%
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
          mutate(NPMSM = map(wts_data,
            ~ survNPMSM(wts_data = .x,
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
          mutate(MSM = map(wts_data,
            ~ survMSM(wts_data = .x,
                      OData = OData,
                      tbreaks = tbreaks,
                      glm_package = "speedglm"))) %>%
          mutate(MSM = map(MSM, "estimates"))

  ## save IPW tables (will be later merged with main results dataset)
  IPWtabs <-  analysis %>%
               left_join(IPW) %>%
               distinct(intervened_TRT, trunc_MSM, wts_tabs, FUPtimes_tabs) %>%
               nest(intervened_TRT, wts_tabs, FUPtimes_tabs, .key = "IPWtabs")

  IPW <- IPW %>% select(-wts_data, -wts_tabs, -FUPtimes_tabs)

  test_that("IPW results w/out Monitoring match", {
    IPW[["NPMSM"]][[1]][["St.NPMSM"]]
     # [1] 0.9896490 0.9896490 0.9896490 0.9896490 0.9896490 0.9747998 0.9747998 0.8488922 0.8488922 0.8488922 0.8488922 0.8488922 0.8340430 0.8340430

    IPW[["NPMSM"]][[2]][["St.NPMSM"]]
    # [1] 0.9896132 0.9794273 0.9794273 0.8814496 0.8182403 0.7827866 0.7378625 0.6762968 0.6762968 0.6318210 0.6318210 0.5800140 0.5800140 0.5800140

    IPW[["MSM"]][[1]][["St.MSM"]]
    # [1] 0.9896490 0.9896490 0.9896490 0.9896490 0.9896490 0.9747998 0.9747998 0.8488922 0.8488922 0.8488922 0.8488922 0.8437729 0.8386845 0.8336267

    IPW[["MSM"]][[1]][["SE.MSM"]]
    # [1] 0.01081701 0.01081701 0.01081701 0.01081701 0.01081701 0.01956379 0.01956379 0.11775100 0.11775100 0.11775100 0.11775100 0.11781627 0.11811526 0.11863771

    IPW[["MSM"]][[2]][["St.MSM"]]
    # [1] 0.9896132 0.9794273 0.9794273 0.8814496 0.8182403 0.7827866 0.7378625 0.6762968 0.6609433 0.6459384 0.6312741 0.6120658 0.5934420 0.5753848

    IPW[["MSM"]][[2]][["SE.MSM"]]
    # [1] 0.01034527 0.01445261 0.01445261 0.07935271 0.08220456 0.08330800 0.08508971 0.08549283 0.08446174 0.08474488 0.08615740 0.08458145 0.08512814 0.08739380
  })
  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  test_that("GCOMP results w/out Monitoring match", {
    GCOMP <-analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule) %>%
          mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
            ~ fit_GCOMP(intervened_TRT = .x,
                          stratifyQ_by_rule = .y,
                          tvals = tvals,
                          OData = OData,
                          # models = models_Q,
                          Qforms = Qforms,
                          fit_method = fit_method_Q
                          ))) %>%
          mutate(GCOMP = map(GCOMP, "estimates"))

    GCOMP[["GCOMP"]][[1]][["St.GCOMP"]]
    # [1] 0.99 0.99 0.99

    GCOMP[["GCOMP"]][[2]][["St.GCOMP"]]
    # [1] 0.9900000 0.9719893 0.9593532
  })

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  test_that("TMLE results w/out Monitoring match", {
    TMLE <- CVTMLE <- analysis %>%
          distinct(intervened_TRT, stratifyQ_by_rule)
    TMLE <- TMLE %>%
          mutate(TMLE = pmap(TMLE, fit_TMLE,
                             tvals = tvals,
                             OData = OData,
                             # models = models_Q,
                             Qforms = Qforms,
                             fit_method = fit_method_Q
                             )) %>%
          mutate(TMLE = map(TMLE, "estimates"))

    TMLE[["TMLE"]][[1]][["St.TMLE"]]
    # [1] 0.99 0.99 0.99
    TMLE[["TMLE"]][[1]][["SE.TMLE"]]
    # [1] 0.009949874 0.009949874 0.009949874

    TMLE[["TMLE"]][[2]][["St.TMLE"]]
    # [1] 0.9900000 0.9789573 0.9679731
    TMLE[["TMLE"]][[2]][["SE.TMLE"]]
    # [1] 0.009949874 0.015098587 0.015938640
  })
