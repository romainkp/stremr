context("Testing sl3 baseline glm learner with various family() args; test sl3 fits with empty strata")
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
  library("tibble")
  library("tidyr")
  library("purrr")
  library("dplyr")

  ## -----------------------------------------------------------------------
  ## Define binomial and gaussian GLM learners
  ## -----------------------------------------------------------------------
  lrn_glm_base_bin <- Lrnr_glm$new(family = binomial())
  ## try new family (as gaussian)
  lrn_glm_base_gaus <- Lrnr_glm$new(family = gaussian())

  ## -----------------------------------------------------------------------
  ## Load data and define some summaries (lags C[t-1], A[t-1], N[t-1])
  ## -----------------------------------------------------------------------
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  # select only the first 100 IDs
  Odat_DT <- Odat_DT[ID %in% (1:500), ]
  setkeyv(Odat_DT, cols = c("ID", "t"))
  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- paste0(lagnodes, ".tminus1")
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

  ## remove indicators of censoring and monitoring events:
  Odat_DT[, "N" := NULL] # Odat_DT[, "C" := NULL]

  ## ------------------------------------------------------------------
  ## Propensity score models
  ## Fit separate model for each strata
  ## ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))

  ## ----------------------s--------------------------------------------------
  ## propensity score (g) w/ PARAMETRIC LOGISTIC REGRESSION (binomial)
  ## ------------------------------------------------------------------------  
  OData <- stremr::importData(Odat_DT, ID = "ID", t_name = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), TRT = "TI", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData,
                         gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT,
                         models_TRT = lrn_glm_base_bin
                        )
  wts_data <- getIPWeights(intervened_TRT = "gTI.dlow", OData = OData)
  mean_gfit_bin <- mean(wts_data[["g0.CAN"]]) # [1] 0.9303566
  mean_IPAW_bin <- mean(wts_data[["cum.IPAW"]]) # [1] 1.085588
  expect_equal(mean_gfit_bin,0.9303566,tolerance = .0000001, scale = 1)
  expect_equal(mean_IPAW_bin,1.085588,tolerance = .000001, scale = 1)

  ## ----------------------s--------------------------------------------------
  ## propensity score (g) w/ PARAMETRIC LINEAR REGRESSION (gaussian)
  ## ------------------------------------------------------------------------  
  OData <- fitPropensity(OData,
                         gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT,
                         models_TRT = lrn_glm_base_gaus
                        )
  wts_data <- getIPWeights(intervened_TRT = "gTI.dlow", OData = OData)
  mean_gfit_gaus <- mean(wts_data[["g0.CAN"]]) # [1] 0.9303652
  mean_IPAW_gaus <- mean(wts_data[["cum.IPAW"]]) # [1] 1.094452
  expect_equal(mean_gfit_gaus,0.9303652,tolerance = .0000001, scale = 1)
  expect_equal(mean_IPAW_gaus,1.094452,tolerance = .000001, scale = 1)

  ## ------------------------------------------------------------------
  ## Test fitting propensity score models with one empty strata
  ## Should do nothing just ignore the strata, same result as before
  ## ------------------------------------------------------------------
  stratify_TRT_nullstrat <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)",                      # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
         "t == 17L"                                            # EMPTY STRATA, FOR TESTING
        ))
  OData_nullstrat <- stremr::importData(Odat_DT, ID = "ID", t_name = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), TRT = "TI", OUTCOME = "Y.tplus1")
  OData_nullstrat <- fitPropensity(OData_nullstrat,
                         gform_TRT = gform_TRT,
                         stratify_TRT = stratify_TRT_nullstrat,
                         models_TRT = lrn_glm_base_bin
                        )
  ## Get the dataset with weights:
  wts_data <- getIPWeights(intervened_TRT = "gTI.dlow", OData = OData_nullstrat)
  mean(wts_data[["g0.CAN"]]) # [1] 0.9303566
  mean(wts_data[["cum.IPAW"]]) # [1] 1.085588
  expect_equal(mean_gfit_bin,mean(wts_data[["g0.CAN"]]))
  expect_equal(mean_IPAW_bin,mean(wts_data[["cum.IPAW"]]))

  