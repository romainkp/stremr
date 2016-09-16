# --------------------------------------------------------------------------------------------------------
# Install h2o (most recent version)
# --------------------------------------------------------------------------------------------------------
# if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
# if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
# # Next, we download packages that H2O depends on.
# if (! ("methods" %in% rownames(installed.packages()))) { install.packages("methods") }
# if (! ("statmod" %in% rownames(installed.packages()))) { install.packages("statmod") }
# if (! ("stats" %in% rownames(installed.packages()))) { install.packages("stats") }
# if (! ("graphics" %in% rownames(installed.packages()))) { install.packages("graphics") }
# if (! ("RCurl" %in% rownames(installed.packages()))) { install.packages("RCurl") }
# if (! ("jsonlite" %in% rownames(installed.packages()))) { install.packages("jsonlite") }
# if (! ("tools" %in% rownames(installed.packages()))) { install.packages("tools") }
# if (! ("utils" %in% rownames(installed.packages()))) { install.packages("utils") }
# # Now we download, install and initialize the H2O package for R.
# install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/rel-turchin/9/R")))
# --------------------------------------------------------------------------------------------------------
# Install data.table (most recent version)
# --------------------------------------------------------------------------------------------------------
# devtools::install_github('Rdatatable/data.table')
# --------------------------------------------------------------------------------------------------------
# Install stremr
# --------------------------------------------------------------------------------------------------------
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)

# ---------------------------------------------------------------------------
# Test h2o randomForest, GBM and deeplearning
# ---------------------------------------------------------------------------
test.h2o.ALL.ML.allestimators10Kdata <- function() {
  `%+%` <- function(a, b) paste0(a, b)
  options(stremr.verbose = TRUE)
  require("data.table")
  require("h2o")

  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
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
  # SET UP h2o cluster & IMPORT DATA
  # ----------------------------------------------------------------
  h2o::h2o.init(nthreads = 1)
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)

  # ----------------------------------------------------------------
  # FIT PROPENSITY SCORES WITH randomForest
  # ----------------------------------------------------------------
  set_all_stremr_options(fit.package = "h2o", fit.algorithm = "randomForest")
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  surv1 <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  surv2 <- survNPMSM(wts.St.dhigh, OData)

  make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                  AddFUPtables = TRUE,
                  openFile = FALSE,
                  WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Jane Doe", y_legend = 0.95)

  # ----------------------------------------------------------------
  # FIT PROPENSITY SCORES WITH gbm
  # ----------------------------------------------------------------
  set_all_stremr_options(fit.package = "h2o", fit.algorithm = "gbm")
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  surv1 <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  surv2 <- survNPMSM(wts.St.dhigh, OData)

  make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                  AddFUPtables = TRUE,
                  openFile = FALSE,
                  WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Jane Doe", y_legend = 0.95)

  # ----------------------------------------------------------------
  # FIT PROPENSITY SCORES WITH deeplearning
  # ----------------------------------------------------------------
  set_all_stremr_options(fit.package = "h2o", fit.algorithm = "deeplearning")
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  surv1 <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  surv2 <- survNPMSM(wts.St.dhigh, OData)

  make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                  AddFUPtables = TRUE,
                  # openFile = FALSE,
                  WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Jane Doe", y_legend = 0.95)

}
