context("NDE assumption")

## ---------------------------------------------------------------------------------
## COMMENTS ON TESTS
## ---------------------------------------------------------------------------------
  require("stremr")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)

  options(width = 100)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
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
  newVarnames <- lagnodes %+% ".tminus1"
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # indicator that the person has never been on treatment up to current t
  Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
  ## ----------------------------------------------------------------
  ## IMPORT DATA
  ## ----------------------------------------------------------------
  set_all_stremr_options(estimator = "speedglm__glm")
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

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

test_that("IPW-KM with stochastic intervention on MONITOR", {
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly")
  surv1.stoch <- survNPMSM(wts.St.dlow, OData)
  expect_true(
  all.equal(
    paste0(surv1.stoch[["estimates"]][["St.NPMSM"]], collapse = ","),
    "0.985853098100471,0.985853098100471,0.985853098100471,0.985853098100471,0.985853098100471,0.978011271542912,0.978011271542912,0.411379895630615,0.411379895630615,0.411379895630615,0.411379895630615,0.411379895630615,0.40998903378949,0.40998903378949,0.409988999668537,0.409988999668537,NaN")
  )

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly")
  surv2.stoch <- survNPMSM(wts.St.dhigh, OData)
  surv2.stoch[["estimates"]][["St.NPMSM"]]
  expect_true(
  all.equal(
    paste0(surv2.stoch[["estimates"]][["St.NPMSM"]], collapse = ","),
    "0.980715607311163,0.974547581455095,0.974547581455095,0.930255027709849,0.926715229154425,0.848554097090501,0.846450767459393,0.734914128610553,0.734914128610553,0.594839291790172,0.594839291790172,0.581772042885852,0.581772042885852,0.581772042885852,0.581772042885852,0.581772042885852,NaN")
  )
})


test_that("IPW-KM with static intervention on MONITOR", {
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101")
  surv1.stat <- survNPMSM(wts.St.dlow, OData)
  expect_true(
  all.equal(
    paste0(surv1.stat[["estimates"]][["St.NPMSM"]], collapse = ","),
    "0.985593308172962,0.985593308172962,0.985593308172962,0.985593308172962,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101")
  surv2.stat <- survNPMSM(wts.St.dhigh, OData)
  expect_true(
  all.equal(
    paste0(surv2.stat[["estimates"]][["St.NPMSM"]], collapse = ","),
    "0.979867574342635,0.94185134182697,0.94185134182697,0.811389084742802,0.811389084742802,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
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
    paste0(surv1.statNDE[["estimates"]][["St.NPMSM"]], collapse = ","),
    "0.98964898762371,0.98964898762371,0.98964898762371,0.98964898762371,0.98964898762371,0.98964898762371,0.98964898762371,0.98964898762371,0.98964898762371,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv2.statNDE <- survNPMSM(wts.St.dhigh, OData)
  expect_true(
  all.equal(
    paste0(surv2.statNDE[["estimates"]][["St.NPMSM"]], collapse = ","),
    "0.989613184270603,0.970013783172476,0.970013783172476,0.936567895363756,0.936567895363756,0.936567895363756,0.936567895363756,0.618726872904148,0.618726872904148,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN")
  )

})


test_that("TMLE / GCOMP with a stochastic intervention on MONITOR", {
  t.surv <- c(0:3)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)

  print(paste0(gcomp_est3[["estimates"]][["St.GCOMP"]], collapse = ","))

  expect_true(
  all.equal(
    paste0(gcomp_est3[["estimates"]][["St.GCOMP"]], collapse = ","),
    "0.989999999994813,0.98351871714058,0.964300896726698,0.901124778851299")
  )

  # stratified modeling by rule followers only:
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = TRUE)

  print(paste0(tmle_est3[["estimates"]][["St.TMLE"]], collapse = ","))

  expect_true(
  all.equal(
    paste0(tmle_est3[["estimates"]][["St.TMLE"]], collapse = ","),
    "0.9900000000139,0.980512185135154,0.971945831548481,0.952551647600147")
  )

  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)

  print(paste0(tmle_est4[["estimates"]][["St.TMLE"]], collapse = ","))

  expect_true(
  all.equal(
    paste0(tmle_est4[["estimates"]][["St.TMLE"]], collapse = ","),
    "0.990000000017604,0.980527554222798,0.973488109369388,0.955062511868786")
  )

})


test_that("TMLE / GCOMP with a static intervention on MONITOR under NDE assumption", {
  t.surv <- c(0:3)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                            useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  expect_true(
  all.equal(
    paste0(gcomp_est3[["estimates"]][["St.GCOMP"]], collapse = ","),
    "0.989999999994757,0.967939608547372,0.961609342301288,0.883919336536042")
  )

  ## stratified modeling by rule followers only:

  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  expect_true(
  all.equal(
    paste0(tmle_est3[["estimates"]][["St.TMLE"]], collapse = ","),
    "0.990000000013414,0.964650728489762,0.943901699485709,0.933132156613137")
  )

  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  expect_true(
  all.equal(
    paste0(tmle_est4[["estimates"]][["St.TMLE"]], collapse = ","),
    "0.990000000017421,0.965155755596925,0.950465241673256,0.869326959247925")
  )

})


test_that("TMLE / GCOMP with a stochastic intervention on MONITOR under NDE assumption", {
  t.surv <- c(0:3)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly",
                            useonly_t_MONITOR = "gPois3.yrly == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  expect_true(
  all.equal(
    paste0(gcomp_est3[["estimates"]][["St.GCOMP"]], collapse = ","),
    "0.989999999994757,0.977605017200815,0.96640066868087,0.900457218199797")
  )

  # stratified modeling by rule followers only:
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly",
                        useonly_t_MONITOR = "gPois3.yrly == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  expect_true(
  all.equal(
    paste0(tmle_est3[["estimates"]][["St.TMLE"]], collapse = ","),
    "0.990000000013414,0.975704302751518,0.964657999663891,0.884770719616591")
  )

  # pooling all observations (no stratification):
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly",
                        useonly_t_MONITOR = "gPois3.yrly == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  expect_true(
  all.equal(
    paste0(tmle_est4[["estimates"]][["St.TMLE"]], collapse = ","),
    "0.990000000017421,0.975597334629591,0.967482584420614,0.887221610893926")
  )

})