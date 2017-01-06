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

