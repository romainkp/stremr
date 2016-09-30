test.directRegressionDefn.10Kdata <- function() {
  options(stremr.verbose = FALSE)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K

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
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)

  # ------------------------------------------------------------------
  # Alternative approach way to specify regression models
  # ------------------------------------------------------------------
  # options(stremr.verbose = TRUE)
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  reg_CENS <- define_single_regression(OData, "C ~ highA1c + t")
  reg_TRT <- c(
      define_single_regression(OData, "TI ~ CVD + highA1c",
          stratify = list(TI = "t == 0L")),
      define_single_regression(OData, "TI ~ CVD + highA1c",
          stratify = list(TI = "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)")),
      define_single_regression(OData, "TI ~ 1",
          stratify = list(TI = "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)")),
      define_single_regression(OData, "TI ~ 1",
          stratify = list(TI = "(t > 0L) & (barTIm1eq0 == 0L)"))
      )
  reg_MONITOR <- define_single_regression(OData, "N ~ 1")
  OData <- fitPropensity(OData, reg_CENS = reg_CENS, reg_TRT = reg_TRT)

  require("magrittr")
  St.dlow2 <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly") %>%
             survNPMSM(OData)  %$%
             IPW_estimates
  # St.dlow2

  St.dhigh2 <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
              survNPMSM(OData) %$%
              IPW_estimates
  # St.dhigh2
}

notest.savedSL.10Kdata <- function() {
  reqh2o <- requireNamespace("h2o", quietly = TRUE)
  reqSL <- requireNamespace("h2oEnsemble", quietly = TRUE)
  if (reqh2o && reqSL) {
    `%+%` <- function(a, b) paste0(a, b)
    require("data.table")
    # options(stremr.verbose = TRUE)
    set_all_stremr_options(fit.package = "h2o", fit.algorithm = "SuperLearner")
    require("h2o")
    require('h2oEnsemble')
    h2o::h2o.init(nthreads = 1)
    # h2o::h2o.init(nthreads = -1)
    # h2o::h2o.init(nthreads = -1, startH2O = FALSE)

    data(OdatDT_10K)
    Odat_DT <- OdatDT_10K

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
    # IMPORT DATA / INIT h2o
    # ----------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)

    # ------------------------------------------------------------------
    # Define models
    # ------------------------------------------------------------------
    glm_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
                             alpha = c(0,1,seq(0.1,0.9,0.1)), lambda = c(0,1e-7,1e-5,1e-3,1e-1))
    h2o.glm.2 <- function(..., x = "highA1c", alpha = 0.0) h2o.glm.wrapper(..., x = x, alpha = alpha)
    learner <- c("h2o.glm.2")
    SLparams = list(fit.package = "h2o", fit.algorithm = "SuperLearner", grid.algorithm = c("glm"), learner = learner, nfolds = 3, seed = 23, glm = glm_hyper_params)
    params_CENS = c(SLparams)
    params_TRT = c(SLparams)
    params_MONITOR = list(fit.package = "speedglm", fit.algorithm = "glm")

    # ------------------------------------------------------------------
    # Fit propensity scores for Treatment, Censoring & Monitoring
    # ------------------------------------------------------------------
    gform_TRT <- c("TI ~ CVD + highA1c")
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
                            params_CENS = params_CENS, params_TRT = params_TRT, params_MONITOR = params_MONITOR)
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
    # Alternative approach way to specify regression models
    # ------------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    reg_CENS <- define_single_regression(OData, "C ~ highA1c + t", params = c(SLparams, save.ensemble = TRUE, ensemble.dir.path = "./h2o-ensemble-model-CENS"))
    reg_TRT <- c(
        define_single_regression(OData, "TI ~ CVD + highA1c",
            stratify = list(TI = "t == 0L"),
            params = c(SLparams, save.ensemble = TRUE, ensemble.dir.path = "./h2o-ensemble-model-TRT1")),
        define_single_regression(OData, "TI ~ CVD + highA1c",
            stratify = list(TI = "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)"),
            params = c(SLparams, save.ensemble = TRUE, ensemble.dir.path = "./h2o-ensemble-model-TRT2")),
        define_single_regression(OData, "TI ~ 1",
            stratify = list(TI = "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)"),
            params = params_MONITOR),
        define_single_regression(OData, "TI ~ 1",
            stratify = list(TI = "(t > 0L) & (barTIm1eq0 == 0L)"),
            params = params_MONITOR)
        )
    reg_MONITOR <- define_single_regression(OData, "N ~ 1", params = params_MONITOR)
    OData <- fitPropensity(OData, reg_CENS = reg_CENS, reg_TRT = reg_TRT, reg_MONITOR = reg_MONITOR)

    require("magrittr")
    St.dlow2 <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly") %>%
               survNPMSM(OData)  %$%
               IPW_estimates
    St.dlow2

    St.dhigh2 <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
                survNPMSM(OData) %$%
                IPW_estimates
    St.dhigh2

    # report.path <- "/Users/olegsofrygin/Dropbox/KP/monitoring_simstudy/stremr_examples"
    # make_report_rmd(OData,
    #                 # MSM = MSM.IPAW,
    #                 # AddFUPtables = TRUE,
    #                 # MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = FALSE),
    #                 # WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
    #                 file.name = "sim.data.example.fup", file.path = report.path, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)

    # ------------------------------------------------------------------
    # USE THE PREVIOUSLY SAVED SL FITS INSTEAD OF FITTIG NEW SL
    # ------------------------------------------------------------------
    OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
    reg_CENS <- define_single_regression(OData, "C ~ highA1c + t", params = c(SLparams, load.ensemble = TRUE, ensemble.dir.path = "./h2o-ensemble-model-CENS"))
    reg_TRT <- c(
        define_single_regression(OData, "TI ~ CVD + highA1c",
            stratify = list(TI = "t == 0L"),
            params = c(SLparams, load.ensemble = TRUE, ensemble.dir.path = "./h2o-ensemble-model-TRT1")),
        define_single_regression(OData, "TI ~ CVD + highA1c",
            stratify = list(TI = "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)"),
            params = c(SLparams, load.ensemble = TRUE, ensemble.dir.path = "./h2o-ensemble-model-TRT2")),
        define_single_regression(OData, "TI ~ 1",
            stratify = list(TI = "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)"),
            params = params_MONITOR),
        define_single_regression(OData, "TI ~ 1",
            stratify = list(TI = "(t > 0L) & (barTIm1eq0 == 0L)"),
            params = params_MONITOR)
        )
    reg_MONITOR <- define_single_regression(OData, "N ~ 1", params = params_MONITOR)

    OData <- fitPropensity(OData, reg_CENS = reg_CENS, reg_TRT = reg_TRT, reg_MONITOR = reg_MONITOR)

    require("magrittr")
    St.dlow3 <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly") %>%
               survNPMSM(OData)  %$%
               IPW_estimates
    St.dlow3

    St.dhigh3 <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
                survNPMSM(OData) %$%
                IPW_estimates
    St.dhigh3

    report.path <- "/Users/olegsofrygin/Dropbox/KP/monitoring_simstudy/stremr_examples"
    make_report_rmd(OData, openFile = FALSE,
                    # MSM = MSM.IPAW,
                    # AddFUPtables = TRUE,
                    # MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = FALSE),
                    # WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                    file.name = "sim.data.example.fup", file.path = report.path, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)

    h2o::h2o.shutdown(prompt = FALSE)
  }
}