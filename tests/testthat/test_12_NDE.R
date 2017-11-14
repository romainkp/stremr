context("NDE assumption")
  run_test <- TRUE

  library("stremr")
  library("magrittr")
  library("data.table")
  library("testthat")
  library("sl3")

  data.table::setDTthreads(1)
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  options(sl3.verbose = FALSE)
  options(condensier.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = TRUE)
  # options(sl3.verbose = TRUE)
  # options(condensier.verbose = TRUE)

  options(width = 100)

  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  Odat_DT <- Odat_DT[ID %in% (1:100), ]
  # define intervention on N as 0101010101...
  Odat_DT[, ("N.star.0101") := t%%2]
  setkeyv(Odat_DT, cols = c("ID", "t"))

  ## ---------------------------------------------------------------------------
  ## Define some summaries (lags C[t-1], A[t-1], N[t-1])
  ## ---------------------------------------------------------------------------
  ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
  lagnodes <- c("C", "TI", "N")
  newVarnames <- paste0(lagnodes, ".tminus1")
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
  ## ----------------------------------------------------------------
  ## IMPORT DATA
  ## ----------------------------------------------------------------
  # set_all_stremr_options(estimator = "speedglm__glm")
  OData <- importData(Odat_DT,
                      ID = "ID", t = "t",
                      covars = c("highA1c", "lastNat1", "lastNat1.factor"),
                      CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome,
                      remove_extra_rows = FALSE)

  ## ------------------------------------------------------------------
  ## Fit propensity scores for Treatment, Censoring & Monitoring
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

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

test_that("IPW-KM with stochastic intervention on MONITOR", {
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly")
  surv1.stoch <- survNPMSM(wts.St.dlow, OData)
  expect_equal(
    surv1.stoch[["estimates"]][["St.NPMSM"]],
    c(0.9858531, 0.9858531, 0.9858531, 0.9858531, 0.9858531, 0.9780113, 0.9780113, 0.4113799, 0.4113799, 0.4113799, 0.4113799, 0.4113799, 0.4099890, 0.4099890, 0.4099890, 0.4099890, NaN),
    tolerance = .0001
  )

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly")
  surv2.stoch <- survNPMSM(wts.St.dhigh, OData)
  expect_equal(
    surv2.stoch[["estimates"]][["St.NPMSM"]],
    c(0.9807156, 0.9745476, 0.9745476, 0.9302550, 0.9267152, 0.8485541, 0.8464508, 0.7349141, 0.7349141, 0.5948393, 0.5948393, 0.5817720, 0.5817720, 0.5817720, 0.5817720, 0.5817720, NaN),
    tolerance = .0001
  )
})

test_that("IPW-KM with static intervention on MONITOR", {
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101")
  surv1.stat <- survNPMSM(wts.St.dlow, OData)
  expect_true(
  all.equal(
    paste0(round(surv1.stat[["estimates"]][["St.NPMSM"]],4), collapse = ","),
    "0.9856,0.9856,0.9856,0.9856,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101")
  surv2.stat <- survNPMSM(wts.St.dhigh, OData)
  expect_true(
  all.equal(
    paste0(round(surv2.stat[["estimates"]][["St.NPMSM"]],4), collapse = ","),
    "0.9799,0.9419,0.9419,0.8114,0.8114,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )
})

test_that("IPW-KM with static intervention on MONITOR under NDE assumption", {
  ## ---------------------------------------------------------------------------------------------------------
  ## IPW-KM with static intervention on MONITOR under NDE assumption
  ## -- Will not intervene on g nodes that satisfy the logical expression provided to this arg
  ## ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv1.statNDE <- survNPMSM(wts.St.dlow, OData)
  expect_true(
  all.equal(
    paste0(round(surv1.statNDE[["estimates"]][["St.NPMSM"]],4), collapse = ","),
    "0.9896,0.9896,0.9896,0.9896,0.9896,0.9896,0.9896,0.9896,0.9896,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv2.statNDE <- survNPMSM(wts.St.dhigh, OData)
  expect_true(
  all.equal(
    paste0(round(surv2.statNDE[["estimates"]][["St.NPMSM"]],4), collapse = ","),
    "0.9896,0.97,0.97,0.9366,0.9366,0.9366,0.9366,0.6187,0.6187,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )
})

test_that("TMLE / GCOMP with a stochastic intervention on MONITOR", {
  t.surv <- c(0:10)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  # print(paste0(gcomp_est3[["estimates"]][["St.GCOMP"]], collapse = ","))
  expect_equal(
    gcomp_est3[["estimates"]][["St.GCOMP"]],
    c(0.9900000, 0.9835187, 0.9643009, 0.9011248, 0.8622261, 0.7791177, 0.7377508, 0.5964215, 0.5882335, 0.5921651, 0.5921651),
    tolerance = .0001
  )

  ## stratified modeling by rule followers only:
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  # print(paste0(tmle_est3[["estimates"]][["St.TMLE"]], collapse = ","))
  ## was failing on bug fixed stremr version:
  # expect_true(
  # all.equal(
  #   paste0(round(tmle_est3[["estimates"]][["St.TMLE"]],4), collapse = ","),
  #   "0.99,0.9805,0.9719,0.9526")
  # )
  expect_equal(
    tmle_est3[["estimates"]][["St.TMLE"]],
    c(0.9900000, 0.9815155, 0.9900000, 0.9492275, 0.9465099, 0.9465099, 0.8165468, 0.6305506, 0.6305506, 0.5363613, 0.5363540),
    tolerance = .0001
  )
  expect_equal(
    tmle_est3[["estimates"]][["SE.TMLE"]],
    c(0.009949874, 0.014289238, 0.013203304, 0.049681300, 0.053829874, 0.097355823, 0.109072012, 0.053244979, 0.053244974, 0.098164546, 0.098165078),
    tolerance = .0001
  )

  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  # print(paste0(tmle_est4[["estimates"]][["St.TMLE"]], collapse = ","))
  ## was failing on bug fixed vs of stremr
  # expect_true(
  # all.equal(
  #   paste0(round(tmle_est4[["estimates"]][["St.TMLE"]],4), collapse = ","),
  #   "0.99,0.9805,0.9735,0.9551")
  # )
  expect_equal(
    tmle_est4[["estimates"]][["St.TMLE"]],
    c(0.9900000, 0.9817518, 0.9744054, 0.9402445, 0.9385821, 0.8250578, 0.8179511, 0.6204416, 0.7057676, 0.6130236, 0.6130236),
    tolerance = .0001
  )
  expect_equal(
    tmle_est4[["estimates"]][["SE.TMLE"]],
    c(0.009949874, 0.014172146, 0.018310629, 0.042300721, 0.046925772, 0.076093885, 0.076334074, 0.101056422, 0.150975762, 0.152185581, 0.152185581),
    tolerance = .0001
  )
})

test_that("TMLE / GCOMP with a static intervention on MONITOR under NDE assumption", {
  t.surv <- c(0:10)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                            useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  expect_equal(
    gcomp_est3[["estimates"]][["St.GCOMP"]],
    c(0.9900000, 0.9679396, 0.9616093, 0.8839193, 0.8246476, 0.8006976, 0.7883042, 0.5805688, 0.4600393, 0.5749507, 0.5749507),
    tolerance = .0001
  )
  ## stratified modeling by rule followers only:
  models_Q <- defModel(estimator = "speedglm__glm")
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE,
                        models = models_Q)

  expect_equal(
    tmle_est3[["estimates"]][["St.TMLE"]],
    c(0.9900000, 0.9646507, 0.9439017, 0.9331322, 0.8400121, 0.8828130, 0.8736917, 0.6038734, 0.6038734, 0.5360617, 0.5360617),
    tolerance = .0001
  )
  expect_equal(
    tmle_est3[["estimates"]][["SE.TMLE"]],
    c(0.009949874, 0.023290794, 0.025426711, 0.049717703, 0.070761799, 0.057952652, 0.056653507, 0.190193465, 0.190193465, 0.197448696, 0.197448696),
    tolerance = .0001
  )
  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  expect_equal(
    tmle_est4[["estimates"]][["St.TMLE"]],
    c(0.9900000, 0.9651558, 0.9504652, 0.8693270, 0.8209692, 0.8604199, 0.8519255, 0.5933615, 0.5262126, 0.5286196, 0.5286196),
    tolerance = .0001
  )
  expect_equal(
    tmle_est4[["estimates"]][["SE.TMLE"]],
    c(0.009949874, 0.023247797, 0.025311490, 0.054319016, 0.067307607, 0.053147718, 0.053000678, 0.173370912, 0.222022409, 0.191511407, 0.191511407),
    tolerance = .0001
  )
})
