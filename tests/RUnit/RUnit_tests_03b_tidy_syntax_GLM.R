test.tidy.speedglm.10Kdata <- function() {
  ## -----------------------------------------------------------------------
  ## ****************************** IMPORTANT ******************************
  ## -----------------------------------------------------------------------
  ## Make sure to always install the latest versions of these packages: "gridisl" and "stremr":
  # devtools::install_github('osofr/gridisl')
  # devtools::install_github('osofr/stremr')
  ## -----------------------------------------------------------------------

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
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = TRUE)
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  set_all_stremr_options(estimator = "speedglm__glm")

  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # select only the first 1,000 IDs
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
  tvals <- 0:3
  tmax <- 13
  tbreaks = c(1:8,11,14)-1

  ## This dataset defines all parameters that we like to vary in this analysis (including different interventions)
  ## That is, each row of this dataset corresponds with a single analysis, for one intervention of interest.
  analysis <- list(intervened_TRT = c("gTI.dlow", "gTI.dhigh"),
                  trunc_wt = c(FALSE, TRUE),
                  stratifyQ_by_rule = c(TRUE, FALSE)) %>%
                  cross_d() %>%
                  arrange(stratifyQ_by_rule) %>%
                  mutate(trunc_MSM = map_dbl(trunc_wt, ~ ifelse(.x, trunc_IPW, Inf))) %>%
                  mutate(trunc_TMLE = trunc_MSM*10)

  ## ----------------------------------------------------------------
  ## IMPORT DATA
  ## ----------------------------------------------------------------
  OData <- stremr::importData(Odat_DT,
                              ID = "ID", t = "t",
                              covars = c("highA1c", "lastNat1", "lastNat1.factor"),
                              CENS = "C", TRT = "TI", MONITOR = "N",
                              OUTCOME = outcome,
                              remove_extra_rows = FALSE)
  OData$dat.sVar[]

  ## ------------------------------------------------------------------------
  ## Define models for iterative G-COMP (Q) -- PARAMETRIC LOGISTIC REGRESSION
  ## ------------------------------------------------------------------------
  ## regression formulas for Q's:
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(tvals)+1))

  ## ----------------------------------------------------------------
  ## Fit propensity score models.
  ## We are using the same model ensemble defined in models_g for censoring, treatment and monitoring mechanisms.
  ## ----------------------------------------------------------------
  OData <- fitPropensity(OData,
                         gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  ## ------------------------------------------------------------
  ## RUN IPW ANALYSES
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

        dplyr::mutate(directIPW = map2(wts_data, trunc_weight,
          ~ directIPW(wts_data = .x,
                      trunc_weights = .y,
                      OData = OData))) %>%
        dplyr::mutate(directIPW = purrr::map(directIPW, "estimates")) %>%
        rename(trunc_MSM = trunc_weight)

  ## save IPW tables (will be later merged with main results dataset)
  IPWtabs <-  analysis %>%
               left_join(IPW) %>%
               distinct(intervened_TRT, trunc_MSM, wts_tabs, FUPtimes_tabs) %>%
               nest(intervened_TRT, wts_tabs, FUPtimes_tabs, .key = "IPWtabs")

  IPW <- IPW %>% select(-wts_data, -wts_tabs, -FUPtimes_tabs)

  IPW[["directIPW"]][[1]]
  attributes(IPW[["directIPW"]][[1]])

  ## ------------------------------------------------------------
  ## GCOMP ANALYSIS
  ## ------------------------------------------------------------
  GCOMP <-analysis %>%
        distinct(intervened_TRT, stratifyQ_by_rule) %>%
        mutate(GCOMP = map2(intervened_TRT, stratifyQ_by_rule,
          ~ fit_GCOMP(intervened_TRT = .x,
                        stratifyQ_by_rule = .y,
                        tvals = tvals,
                        OData = OData,
                        Qforms = Qforms))) %>%
        mutate(GCOMP = map(GCOMP, "estimates"))

  ## ------------------------------------------------------------
  ## TMLE ANALYSIS
  ## ------------------------------------------------------------
  TMLE <- CVTMLE <- analysis %>%
          rename(trunc_weight = trunc_TMLE) %>%
          distinct(intervened_TRT, stratifyQ_by_rule, trunc_weight)

  TMLE <- TMLE %>%
        mutate(TMLE = pmap(TMLE, fit_TMLE,
                           tvals = tvals,
                           OData = OData,
                           Qforms = Qforms)) %>%
        mutate(TMLE = map(TMLE, "estimates")) %>%
        rename(trunc_TMLE = trunc_weight)

  # CVTMLE <- CVTMLE %>%
  #       mutate(CVTMLE = pmap(CVTMLE, fit_CVTMLE,
  #                          tvals = tvals,
  #                          OData = OData,
  #                          Qforms = Qforms)) %>%
  #       mutate(CVTMLE = map(CVTMLE, "estimates")) %>%
  #       rename(trunc_TMLE = trunc_weight)

  ## ------------------------------------------------------------
  ## COMBINE ALL ANALYSES INTO A SINGLE DATASET
  ## ------------------------------------------------------------
  results <-  analysis %>%
              left_join(IPW) %>%
              left_join(GCOMP) %>%
              left_join(TMLE)
              #  %>%
              # left_join(CVTMLE)

  ## Nest each estimator by treatment regimen (we now only show the main analysis rows)
  results <- results %>%
              # , CVTMLE
             nest(intervened_TRT, NPMSM, MSM.crude, MSM, directIPW, GCOMP, TMLE, .key = "estimates")

  ## Calculate RDs (contrasting all interventions, for each analysis row & estimator).
  ## The RDs data no longer needs the intervened_TRT column
  results <-  results %>%
              mutate(RDs =
                map(estimates,
                  ~ select(.x, -intervened_TRT) %>%
                  map(~ get_RDs(.x)) %>%
                  as_tibble()
                  ))

}