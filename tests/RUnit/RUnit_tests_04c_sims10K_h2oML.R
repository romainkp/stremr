# --------------------------------------------------------------------------------------------------------
# Install data.table (most recent version)
# devtools::install_github('Rdatatable/data.table')
# --------------------------------------------------------------------------------------------------------
## For installing most recent vs of h2o see: https://s3.amazonaws.com/h2o-release/h2o/master/latest.html
# --------------------------------------------------------------------------------------------------------
# Install stremr:
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test h2o randomForest, GBM and deeplearning
# ---------------------------------------------------------------------------
test.h2o.ALL.ML.allestimators10Kdata <- function() {
    reqh2o <- requireNamespace("h2o", quietly = TRUE)
    if (!reqh2o) return(NULL)

    options(width = 100)
    `%+%` <- function(a, b) paste0(a, b)
    require("data.table")
    require("h2o")
    require("sl3")
    library("stremr")
    # library("gridisl")
    # options(stremr.verbose = TRUE)
    # options(sl3.verbose = TRUE)
    # options(gridisl.verbose = TRUE)
    options(stremr.verbose = FALSE)
    options(gridisl.verbose = FALSE)

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K
    # select only the first 100 IDs
    Odat_DT <- Odat_DT[ID %in% (1:100), ]
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
    # FIT PROPENSITY SCORES WITH gbm
    # ----------------------------------------------------------------
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR, 
                            models_TRT = sl3::Lrnr_h2o_grid$new(algorithm = "gbm", ntrees = 10))

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv1 <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv2 <- survNPMSM(wts.St.dhigh, OData)

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    openFile = FALSE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ----------------------------------------------------------------
    # FIT PROPENSITY SCORES WITH deeplearning
    # ----------------------------------------------------------------
    OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                            stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR, 
                            models_TRT = sl3::Lrnr_h2o_grid$new(algorithm = "deeplearning"))

    wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
    surv1 <- survNPMSM(wts.St.dlow, OData)
    wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
    surv2 <- survNPMSM(wts.St.dhigh, OData)

    if (rmarkdown::pandoc_available(version = "1.12.3"))
        make_report_rmd(OData, NPMSM = list(surv1, surv2), wts_data = list(wts.St.dlow, wts.St.dhigh),
                    AddFUPtables = TRUE,
                    openFile = FALSE,
                    WTtables = get_wtsummary(list(wts.St.dlow, wts.St.dhigh), cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ h2o random forest, using gridisl
    # ---------------------------------------------------------------------------------------------------------
    params = gridisl::defModel(estimator = "h2o__randomForest", ntrees = 10, learn_rate = 0.1, sample_rate = 0.9, col_sample_rate = 0.9, balance_classes = TRUE)
    # params = list(fit.package = "h2o", fit.algorithm = "randomForest", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)
    t.surv <- c(2)
    Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)
    tmle_est$estimates[]

    # ---------------------------------------------------------------------------------------------------------
    # TMLE w/ h2o random forest, using sl3 (should be equivalent)
    # ---------------------------------------------------------------------------------------------------------
    # params = sl3::Lrnr_h2o_grid$new(algorithm = "randomForest", ntrees = 10, learn_rate = 0.1, sample_rate = 0.9, col_sample_rate = 0.9, balance_classes = TRUE)
    # # params = list(fit.package = "h2o", fit.algorithm = "randomForest", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)
    # t.surv <- c(2)
    # Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
    # tmle_est_sl3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)
    # tmle_est_sl3$estimates[]
    # checkEquals(tmle_est$estimates[], tmle_est_sl3$estimates[])
    
    h2o::h2o.shutdown(prompt = FALSE)

}