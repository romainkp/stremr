  require("data.table")
  require("stremr")
  options(stremr.verbose = TRUE)
  options(gridisl.verbose = FALSE)

  # ------------------------------------------------------------------------------------------------------
  # load simulated data
  # ------------------------------------------------------------------------------------------------------
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # Odat_DT <- Odat_DT[ID %in% (1:100), ]
  setkeyv(Odat_DT, cols = c("ID", "t"))

  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- lagnodes %+% ".tminus1"
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
                              # , MONITOR = "N",
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
  ## this is your \bar{a} identifier of the static regimen 1, i.e., f(\bar{a}) create your own variable
  wts.DT.0[, "theta" := 0]
  wts.DT.0[, "sumtheta" := cumsum(theta), by = ID]

  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "gTI.dhigh") # , rule_name = "TI0"
  ## this is your \bar{a} identifier of the static regimen 2, i.e., f(\bar{a}) create your own variable
  wts.DT.1[, "theta" := 1]
  wts.DT.1[, "sumtheta" := cumsum(theta), by = ID]

  # ----------------------------------------------------------------------
  # old IPW-MSM for hazard, only allows smoothing over time-intervals (fit constant hazard on some range of t values)
  # ----------------------------------------------------------------------
  ## old MSM that allows pooling over time intervals, summaries f(t,d)=(theta, sumtheta) aren't allowed
  survMSM_res_old <- survMSM(list(wts.DT.1, wts.DT.0), OData, glm_package = "speedglm", return_wts = TRUE)
  survMSM_res_old[["gTI.dlow"]]
  survMSM_res_old[["gTI.dhigh"]]

  # ----------------------------------------------------------------------
  ## new IPW-MSM for the hazard allows using arbitrary GLM formulas
  # ----------------------------------------------------------------------
  ## EXAMPLE 1.
  ## Saturated MSM, must match the MSM above (define a separate indicator for each (t x rule) combo)
  survMSM_res <- fit_hMSM(list(wts.DT.1, wts.DT.0),
                           form = "Y.tplus1 ~ -1 + as.factor(t):as.factor(rule.name)",
                           OData)
  survMSM_res[["estimates"]]
  survMSM_res[["msm.fit"]]

  ## EXAMPLE 2.
  ## MSM with user-defined covariate f(t,d)='theta' that is contant in t (only varies by d)
  ## This MSM will smooth over time (time used as a continuous covariate) and adds interaction of (t,theta)
  survMSM_res2 <- fit_hMSM(list(wts.DT.1, wts.DT.0), form = "Y.tplus1 ~ -1 + t + theta + t:theta",
                            OData,
                            tmax = 10)
  survMSM_res2[["estimates"]]
  survMSM_res2[["msm.fit"]]

  ## EXAMPLE 3.
  ## MSM with user-defined covariate f(t,d)='sumtheta' that varies in (t,d), along with interaction (theta,sumtheta)
  ## add an arbitrary covariate that var
  survMSM_res3 <- fit_hMSM(list(wts.DT.1, wts.DT.0), form = "Y.tplus1 ~ -1 + t + theta + sumtheta", OData)
  survMSM_res3[["estimates"]]
  survMSM_res3[["msm.fit"]]

