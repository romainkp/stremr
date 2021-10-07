context("IPW-MSM for hazard")
  require("data.table")
  require("stremr")
  library("sl3")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  options(sl3.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = TRUE)
  # options(sl3.verbose = TRUE)

  # ------------------------------------------------------------------------------------------------------
  # load simulated data
  # ------------------------------------------------------------------------------------------------------
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  Odat_DT <- Odat_DT[ID %in% (1:500), ]
  setkeyv(Odat_DT, cols = c("ID", "t"))
  Odat_DT <- copy(Odat_DT)

  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- paste0(lagnodes, ".tminus1")
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

  ## ------------------------------------------------------------------
  ## propensity score models for Treatment, Censoring & Monitoring
  ## ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))

  gform_CENS <- c("C ~ highA1c + t")

  OData <- stremr::importData(Odat_DT,
                              ID = "ID", t = "t",
                              covars = c("highA1c", "lastNat1", "lastNat1.factor"),
                              CENS = "C", TRT = "TI",
                              OUTCOME = outcome,
                              remove_extra_rows = FALSE)

  OData <- fitPropensity(OData,
                         gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT) # , gform_MONITOR = gform_MONITOR

  # ----------------------------------------------------------------------
  # Evaluate the IPW weights for each regimen,
  ## define new summaries f(d,t) = ("theta", "sumtheta") for the hazard MSM
  ## - d is the unique regimen identifier
  ## - t is the time identifier
  ## Note that time variable and "rule.name" will be always available by default.
  # ----------------------------------------------------------------------
  wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "gTI.dlow") # , rule_name = "TI1"
  wts.DT.0 <- copy(wts.DT.0)
  ## this is your \bar{a} identifier of the static regimen 1, i.e., f(\bar{a}) create your own variable
  wts.DT.0[, "theta" := 0]
  wts.DT.0[, "sumtheta" := cumsum(theta), by = ID]

  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "gTI.dhigh") # , rule_name = "TI0"
  wts.DT.1 <- copy(wts.DT.1)
  ## this is your \bar{a} identifier of the static regimen 2, i.e., f(\bar{a}) create your own variable
  wts.DT.1[, "theta" := 1]
  wts.DT.1[, "sumtheta" := cumsum(theta), by = ID]


  # ----------------------------------------------------------------------
  # old IPW-MSM for hazard, only allows smoothing over time-intervals (fit constant hazard on some range of t values)
  # ----------------------------------------------------------------------
  test_that(" old MSM that only allows pooling over time intervals", {
    ## old MSM that allows pooling over time intervals, summaries f(t,d)=(theta, sumtheta) aren't allowed
    survMSM_res_old <- survMSM(list(wts.DT.1, wts.DT.0), OData, glm_package = "speedglm", return_wts = TRUE)
    survMSM_res_old[["gTI.dlow"]][["estimates"]]
    survMSM_res_old[["gTI.dhigh"]][["estimates"]]
  })


  # ----------------------------------------------------------------------
  ## new IPW-MSM for the hazard allows using arbitrary GLM formulas
  # ----------------------------------------------------------------------
  test_that("EXAMPLE 1. saturated hazard MSM, must match the MSM above (define a separate indicator for each (t x rule) combo)", {
  ## EXAMPLE 1.
  ## Saturated MSM, must match the MSM above (define a separate indicator for each (t x rule) combo)
  survMSM_res <- fit_hMSM(list(wts.DT.1, wts.DT.0),
                           form = "Y.tplus1 ~ -1 + as.factor(t):as.factor(rule.name)",
                           OData)
  survMSM_res[["estimates"]]
  survMSM_res[["msm.fit"]]
  })


  test_that("EXAMPLE 2. hazard MSM with user-defined covariate f(t,d)='theta' that is contant in t (only varies by d)", {
    ## EXAMPLE 2.
    ## MSM with user-defined covariate f(t,d)='theta' that is contant in t (only varies by d)
    ## This MSM will smooth over time (time used as a continuous covariate) and adds interaction of (t,theta)
    survMSM_res2 <- fit_hMSM(list(wts.DT.1, wts.DT.0), form = "Y.tplus1 ~ -1 + t + theta + t:theta",
                              OData,
                              tmax = 10)
    survMSM_res2[["estimates"]]
    survMSM_res2[["msm.fit"]]
  })


  test_that("EXAMPLE 3. hazard MSM with user-defined covariate f(t,d)='sumtheta' that varies in (t,d), along with interaction (theta,sumtheta)", {
    ## EXAMPLE 3.
    ## MSM with user-defined covariate f(t,d)='sumtheta' that varies in (t,d), along with interaction (theta,sumtheta)
    ## add an arbitrary covariate that var
    survMSM_res3 <- fit_hMSM(list(wts.DT.1, wts.DT.0), form = "Y.tplus1 ~ -1 + t + theta + sumtheta", OData)
    survMSM_res3[["estimates"]]
    survMSM_res3[["msm.fit"]]
  })


  test_that("EXAMPLE 4. new (non-hazard) MSM with user-defined covariate f(t,d)='sumtheta' that varies in (t,d), along with interaction (theta,sumtheta)", {
    ## EXAMPLE 3.
    ## MSM with user-defined covariate f(t,d)='sumtheta' that varies in (t,d), along with interaction (theta,sumtheta)
    ## can add any arbitrary bsl covariate V
    MSM_res4 <- fit_gMSM(list(wts.DT.1, wts.DT.0),
                         ## Adding V to MSM formula:
                         form = Y.tplus1 ~ t + I(rule.name=='gTI.dhigh') + CVD + highA1c,
                         family = "quasibinomial",
                         OData = OData,
                         use_weights = TRUE,
                         stabilize = FALSE,
                         tvals = c(10, 11, 12, 13, 14),
                         # tmax = 12,
                         verbose = TRUE)

    MSM_res4[["msm.fit"]]
    MSM_res4[["beta.SE"]]

   # form = "Y.tplus1 ~ -1 + t + theta + sumtheta", 
   # form = "Y.tplus1 ~ t + theta + sumtheta",
   # form = "Y.tplus1 ~ -1 + t + theta + t:theta",
   # form = "Y.tplus1 ~ -1 + as.factor(t):as.factor(rule.name)"
   # form = "Y.tplus1 ~ t + I(rule.name=='gTI.dhigh')",
  })


# (simpleMSM <- shifted.OUTCOME%+%"~ intnum"%+%" + "%+%paste(paste0("I(rule.name=='",rule.names[rule.names!=dRef],"')"),collapse=" + "))
# HR.MSM.KM <- fit_hMSM(wts_data=IPW.by.rule, form=simpleMSM, OData= OData, use_weights = FALSE, stabilize = TRUE, trunc_weights=Inf, getSEs = TRUE, return_wts = FALSE, tmax = tmax, verbose = FALSE)
# HR.stabIPW.MSM.KM <- fit_hMSM(wts_data=IPW.by.rule, form=simpleMSM, OData= OData, use_weights = TRUE, stabilize = TRUE, trunc_weights=Inf, getSEs = TRUE, return_wts = FALSE, tmax = tmax, verbose = TRUE)



