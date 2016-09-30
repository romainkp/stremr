test.GCOMP.TMLE.10Kdata <- function() {
  options(stremr.verbose = FALSE)
  `%+%` <- function(a, b) paste0(a, b)
  # ---------------------------------------------------------------------------
  # INSTALL CORRECT VERSIONS of data.table and stremr from github:
  # ---------------------------------------------------------------------------
  # devtools::install_github('Rdatatable/data.table')
  # devtools::install_github('osofr/stremr', build_vignettes = FALSE)
  require("data.table")

  # # ---------------------------------------------------------------------------
  # # Test data set included in stremr:
  # # ---------------------------------------------------------------------------
  # head(O.data)
  data(OdatDT_10K)
  ID <- "ID"; t <- "t"; TRT <- "TI"; CENS <- "C"; MONITOR <- "N"; outcome <- "Y.tplus1"; I <- "highA1c";

  # ---------------------------------------------------------------------------
  # DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
  # Might expand this in the future to allow defining arbitrary summaries
  # ---------------------------------------------------------------------------
  # Odat_DT <- obsDTg05_1mil
  Odat_DT <- OdatDT_10K
  lagnodes <- c("C", "TI", "N")
  newVarnames <- lagnodes %+% ".tminus1"
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # Indicator that the person has never been on treatment up to current t
  Odat_DT[, "barTIm1eq0" := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
  # Odat_DT[1:100, ]

  # --------------------------------
  # Define global options for stremr (which R packages to use for model fitting)
  # --------------------------------
  # options(stremr.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")

  # import data into stremr object:
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)

  # --------------------------------
  # Fitting the propensity scores for observed variables (A,C,N)
  # --------------------------------
  # + N.tminus1
  gform_TRT <- "TI ~ CVD + highA1c"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))
  gform_CENS <- c("C ~ highA1c")
  stratify_CENS <- list(C=c("t < 16", "t == 16"))
  gform_MONITOR <- "N ~ 1"

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT,
                                stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  # get IPW-adjusted and KM survival (with hazards over time)
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  St.dlow <- survNPMSM(wts.St.dlow, OData)
  St.dlow

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  St.dhigh <- survNPMSM(wts.St.dhigh, OData)
  St.dhigh

  # ---------------------------------------------------------------------------------------------------------
  # GCOMP AND TMLE w/ GLMs
  # ---------------------------------------------------------------------------------------------------------
  # t.surv <- c(0,1,2,3,4,5,6,7,8,9,10)
  t.surv <- c(1,2,3,10)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  # stratified modeling by rule followers only:
  gcomp_est1 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  names(gcomp_est1)
  gcomp_est1$estimates[]
  gcomp_est1$est_name
  gcomp_est1$IC.Var.S.d
  names(tmle_est1)
  tmle_est1$estimates[]
  tmle_est1$est_name
  tmle_est1$IC.Var.S.d

  # pooling all observations (no stratification):
  gcomp_est2 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est2 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est2$estimates[]; tmle_est2$estimates[]

  # stratified modeling by rule followers only:
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  gcomp_est3$estimates[]; tmle_est3$estimates[]

  # pooling all observations (no stratification):
  gcomp_est4 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est4$estimates[]; tmle_est4$estimates[]

  # ------------------------------------------------------------------------
  # RUN PARALLEL seq-GCOMP & TMLE over t.surv (MUCH FASTER)
  # ------------------------------------------------------------------------
  # require("doParallel")
  # registerDoParallel(cores = 2)
  # data.table::setthreads(1)

  # gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # gcomp_est; tmle_est

  # gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # gcomp_est; tmle_est

  # ---------------------------------------------------------------------------------------------------------
  # GCOMP AND TMLE w/ h2o random forest
  # ---------------------------------------------------------------------------------------------------------
  # require("h2o")
  # h2o::h2o.init(nthreads = 4)
  # # h2o::h2o.init()
  # params = list(fit.package = "h2o", fit.algorithm = "RF", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)

  # t.surv <- c(1,2,3,4,5,6,7,8,9,10)
  # Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  # gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  # tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  # gcomp_est; tmle_est

  # gcomp_fit <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  # tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  # gcomp_est; tmle_est

  # ------------------------------------------------------------------------
  # TEST FOR PROBLEMS WITH > 1 REGRESSION AND > 1 STATA (SHOULD WORK)
  # ------------------------------------------------------------------------
  gform_CENS_test <- c("C1 ~ highA1c", "C2 ~ highA1c")
  stratify_CENS_test <- list(C1=c("t < 16", "t == 16"), C2=c("t < 16", "t == 16"))
  Odat_DT_test <- Odat_DT
  Odat_DT_test[, "C1" := C]
  Odat_DT_test[, "C2" := C]

  OData <- importData(Odat_DT_test, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = c("C1","C2"), TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  OData <- fitPropensity(OData, gform_CENS = gform_CENS_test, stratify_CENS = stratify_CENS_test, gform_TRT = gform_TRT,
                                stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

}