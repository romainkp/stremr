## --------------------------------------------------------------------------------------------------------
## Install data.table (most recent version)
## --------------------------------------------------------------------------------------------------------
# devtools::install_github('Rdatatable/data.table')
## --------------------------------------------------------------------------------------------------------
## For installing most recent vs of h2o see: https://s3.amazonaws.com/h2o-release/h2o/master/latest.html
## --------------------------------------------------------------------------------------------------------
# if ("package:h2o" %in% search()) detach("package:h2o", unload=TRUE)
# if ("h2o" %in% rownames(installed.packages())) remove.packages("h2o")
## Next, download H2O package dependencies:
# pkgs <- c("methods","statmod","stats","graphics","RCurl","jsonlite","tools","utils")
# new.pkgs <- setdiff(pkgs, rownames(installed.packages()))
# if (length(new.pkgs)) install.packages(new.pkgs)
## Download and install the H2O package for R:
# install.packages("h2o", type="source", repos=(c("https://s3.amazonaws.com/h2o-release/h2o/master/3636/R")))
## --------------------------------------------------------------------------------------------------------
## Install h2oEnsemble (most recent stable version 1.8)
## --------------------------------------------------------------------------------------------------------
# install.packages("https://h2o-release.s3.amazonaws.com/h2o-ensemble/R/h2oEnsemble_0.1.8.tar.gz", repos = NULL)
# # Install h2oEnsemble (dev version):
# devtools::install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")
## --------------------------------------------------------------------------------------------------------
## Install stremr
## --------------------------------------------------------------------------------------------------------
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)

# ---------------------------------------------------------------------------
# Test speedglm
# ---------------------------------------------------------------------------
test.speedglm.allestimators10Kdata <- function() {
  options(stremr.verbose = FALSE)
  # options(stremr.verbose = TRUE)
  options(width = 100)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  # obsDTg05_500K <- data.table::fread(input = "./obsDTg05_500K.csv", header = TRUE, na.strings = "NA_h2o")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  setkeyv(Odat_DT, cols = c("ID", "t"))
  head(Odat_DT)
  nrow(Odat_DT)

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

  # ----------------------------------------------------------------
  # IMPORT DATA
  # ----------------------------------------------------------------
  # options(stremr.verbose = TRUE)
  # set_all_stremr_options(fit.package = "glm", fit.algorithm = "glm")
  set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  # To inspect the input data.table:
  # OData$dat.sVar

  # ------------------------------------------------------------------
  # Fit propensity scores for Treatment, Censoring & Monitoring
  # ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))

  gform_CENS <- c("C ~ highA1c + t")
  # stratify_CENS <- list(C=c("t < 16", "t == 16"))
  # stratify_CENS <- list()

  gform_MONITOR <- "N ~ 1"
  # **** really want to define it like this ****
  # gform_TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
  #               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  surv1 <- survNPMSM(wts.St.dlow, OData)

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  surv2 <- survNPMSM(wts.St.dhigh, OData)

  # ------------------------------------------------------------------
  # Piping the workflow
  # ------------------------------------------------------------------
  require("magrittr")
  St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly") %>%
             survNPMSM(OData)  %$%
             IPW_estimates
  St.dlow

  St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
              survNPMSM(OData) %$%
              IPW_estimates
  St.dhigh

  # ------------------------------------------------------------------
  # Running IPW-adjusted MSM for the hazard
  # ------------------------------------------------------------------
  MSM.IPAW <- survMSM(OData,
                      wts_data = list(dlow = wts.St.dlow, dhigh = wts.St.dhigh),
                      t_breaks = c(1:8,12,16)-1,
                      est_name = "IPAW", getSEs = TRUE)
  # names(MSM.IPAW)
  # nrow(t(MSM.IPAW$IC.Var.S.d$gTI.dlow$IC.S))
  # nrow(t(MSM.IPAW$IC.Var.S.d$gTI.dhigh$IC.S))
  # ncol(t(MSM.IPAW$IC.Var.S.d$gTI.dlow$IC.S))
  # ncol(t(MSM.IPAW$IC.Var.S.d$gTI.dhigh$IC.S))
  # # MSM.IPAW$St
  # names(MSM.IPAW$wts_data)
  # OData$nuniqueIDs

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")

  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]

  # ------------------------------------------------------------------------
  # RUN in PARALLEL seq-GCOMP & TMLE over t.surv (MUCH FASTER)
  # ------------------------------------------------------------------------
  # require("doParallel")
  # registerDoParallel(cores = 2)
  if (exists("setthreads")) data.table::setthreads(1)

  t.surv <- c(0,1,4)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est1 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "pooledGCOMP.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  gcomp_est1$estimates[]

  gcomp_est2 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "pooledGCOMP.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  gcomp_est2$estimates[]

  # tmle_est_par1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "pool.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE, parallel = TRUE)
  tmle_est_par1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "pooledTMLE.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE, parallel = FALSE)
  tmle_est_par1$estimates[]
  # tmle_est_par2 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "pool.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE, parallel = TRUE)
  tmle_est_par2 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "pooledTMLE.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE, parallel = FALSE)
  tmle_est_par2$estimates[]

  # tmle_est_par3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "strat.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE, parallel = TRUE)
  tmle_est_par3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "stratTMLE.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE, parallel = FALSE)
  tmle_est_par3$estimates[]
  # tmle_est_par4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "strat.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE, parallel = TRUE)
  tmle_est_par4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "stratTMLE.dlow", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE, parallel = FALSE)
  tmle_est_par4$estimates[]


  # ------------------------------------------------------------------
  # Make a report:
  # ------------------------------------------------------------------
  # report.path <- "/home/ubuntu/stremr_example"
  # file.path = report.path,
  # test for opening file in local OS
  make_report_rmd(OData, file.name = "sim.data.example.fup2", title = "Custom", author = "Jane Doe", openFile = FALSE)

  make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  AddFUPtables = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Author Name", y_legend = 0.99, x_legend = 9.5)

  make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "pdf",
                  AddFUPtables = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Author Name", y_legend = 0.99, x_legend = 9.5)

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a stochastic intervention on MONITOR
  # *** make a separate function
  # **** +1 add evaluation of IPW & TMLE under static NDE
  # **** +2 add evaluation of IPW & TMLE under stochastic NDE intervention on monitoring
  # **** +3 add evaluation of IPW & TMLE under stochastic intervention on treatment
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]
}


# ---------------------------------------------------------------------------
# Test speedglm
# ---------------------------------------------------------------------------
test.speedglm.stochastic.TMLE.NDE.1Kdata <- function() {
  options(stremr.verbose = FALSE)
  options(width = 100)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  Odat_DT <- Odat_DT[ID %in% (1:5000), ]
  # define intervention on N as 0101010101...
  Odat_DT[, ("N.star.0101") := t%%2]
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
  # ----------------------------------------------------------------
  # IMPORT DATA
  # ----------------------------------------------------------------
  set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  # ------------------------------------------------------------------
  # Fit propensity scores for Treatment, Censoring & Monitoring
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

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  # ---------------------------------------------------------------------------------------------------------
  # IPW-KM with stochastic intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly")
  surv1.stoch <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly")
  surv2.stoch <- survNPMSM(wts.St.dhigh, OData)

  # ---------------------------------------------------------------------------------------------------------
  # IPW-KM with static intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101")
  surv1.stat <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101")
  surv2.stat <- survNPMSM(wts.St.dhigh, OData)

  # ---------------------------------------------------------------------------------------------------------
  # IPW-KM with static intervention on MONITOR under NDE assumption
  # ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv1.statNDE <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv2.statNDE <- survNPMSM(wts.St.dhigh, OData)

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a stochastic intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a static intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a static intervention on MONITOR under NDE assumption
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                            useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]
}
