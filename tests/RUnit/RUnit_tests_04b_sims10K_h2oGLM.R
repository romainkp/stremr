# --------------------------------------------------------------------------------------------------------
# Install data.table (most recent version)
# devtools::install_github('Rdatatable/data.table')
# --------------------------------------------------------------------------------------------------------
## For installing most recent vs of h2o see: https://s3.amazonaws.com/h2o-release/h2o/master/latest.html
# --------------------------------------------------------------------------------------------------------
# Install stremr:
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# ---------------------------------------------------------------------------

test.h2oglm.IPW.MSM.10Kdata <- function() {
  reqh2o <- requireNamespace("h2o", quietly = TRUE)
  reqxgb <- requireNamespace("xgboost", quietly = TRUE)
  if (!reqh2o || !reqxgb) { return(NULL) }

    options(width = 100)
    `%+%` <- function(a, b) paste0(a, b)
    require("data.table")
    require("h2o")
    library("stremr")
    options(stremr.verbose = TRUE)
    options(gridisl.verbose = TRUE)
    # options(stremr.verbose = FALSE)
    # options(gridisl.verbose = FALSE)
    # set_all_stremr_options(estimator = "speedglm__glm")

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    Odat_DT <- Odat_DT[ID %in% (1:50), ]
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
    # options(stremr.verbose = TRUE)
    # set_all_stremr_options(estimator = "xgboost__glm", fit_method = "cv", fold_column = "fold_ID")
    # set_all_stremr_options(estimator = "h2o__glm", fit_method = "cv", fold_column = "fold_ID")
    # h2o::h2o.init(nthreads = 1)

    OData <- importData(Odat_DT,
                        ID = "ID", t = "t",
                        covars = c("highA1c", "lastNat1", "lastNat1.factor"),
                        CENS = "C", TRT = "TI", MONITOR = "N",
                        OUTCOME = outcome)
    # to see the input data.table:
    # OData$dat.sVar
    OData <- define_CVfolds(OData, nfolds = 3, fold_column = "fold_ID", seed = 12345)
    OData$dat.sVar[]
    OData$fold_column <- NULL
    OData$nfolds <- NULL

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
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                            family = "quasibinomial",
                            # family = "binomial", solver = "L_BFGS", lambda_search = FALSE,
                            estimator = "xgboost__glm", 
                            fit_method = "cv", fold_column = "fold_ID", nthread = 1)
    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")

    # ------------------------------------------------------------------
    # Running IPW-adjusted MSM for the hazard
    # ------------------------------------------------------------------
    MSM.IPAW <- survMSM(OData,
                        wts_data = list(dlow = wts.St.dlow, dhigh = wts.St.dhigh),
                        tbreaks = c(1:8,12,16)-1,
                        est_name = "IPAW", getSEs = TRUE)

    # names(MSM.IPAW)
    # MSM.IPAW[["St"]]
    # MSM.IPAW[["estimates"]]
    # pl <- ggsurv(MSM.IPAW[["estimates"]])
    # pl

    # names(MSM.IPAW)
    # MSM.IPAW$St
    # if (rmarkdown::pandoc_available(version = "1.12.3"))
        # make_report_rmd(OData,
        #             MSM = MSM.IPAW,
        #             # AddFUPtables = TRUE,
        #             # openFile = FALSE,
        #             # format="pdf",
        #             # openFile = TRUE,
        #             openFile = FALSE,
        #             RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
        #             WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
        #             file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

}
