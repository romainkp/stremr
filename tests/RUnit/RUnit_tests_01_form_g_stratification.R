# --------------------------------
# TEST HELPER FUNCTIONS FOR ID NAMES & calls "DT[,,by=ID]" (frequently does't work when ID="ID")
# --------------------------------
test.helperfuns <- function() {
  require("data.table")
  # O.data <- simulateDATA.fromDAG(Nsize = 1000, rndseed = 124356)
  data(OdataNoCENS)
  # head(OdataNoCENS)

  OdataNoCENS[OdataNoCENS[,"t"]%in%16,"lastNat1"] <- NA
  OdataNoCENS <- OdataNoCENS[,!names(OdataNoCENS)%in%c("highA1c.UN", "timelowA1c.UN")]
  # head(OdataNoCENS)
  OdataNoCENS.DT <- as.data.table(OdataNoCENS, key=c("ID", "t"))
  addN.t1 <- defineMONITORvars(OdataNoCENS, ID = "ID", t = "t", imp.I = "N",
                        MONITOR.name = "N.new", tsinceNis1 = "last.Nt")
  # addN.t1[]
  addN.t2 <- defineMONITORvars(OdataNoCENS.DT, ID = "ID", t = "t", imp.I = "N",
                               MONITOR.name = "N.new", tsinceNis1 = "last.Nt")
  # addN.t2[]
  OdataNoCENS_dhigh_dlow1 <- defineIntervedTRT(OdataNoCENS, theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
                                      CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1", new.TRT.names = c("dlow", "dhigh"))
  # OdataNoCENS_dhigh_dlow1[]

  OdataNoCENS_dhigh_dlow2 <- defineIntervedTRT(OdataNoCENS.DT, theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
                                        CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1", new.TRT.names = c("dlow", "dhigh"))
  # OdataNoCENS_dhigh_dlow2[]

  setnames(OdataNoCENS.DT,old = "ID",new = "ID.expression")
  addN.t1 <- defineMONITORvars(OdataNoCENS.DT, ID = "ID.expression", t = "t", imp.I = "N", MONITOR.name = "N.new", tsinceNis1 = "last.Nt")
  # addN.t1[]

  OdataNoCENS_dhigh_dlow1 <- defineIntervedTRT(OdataNoCENS.DT, theta = c(0,1), ID = "ID.expression", t = "t", I = "highA1c",
                                        CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1", new.TRT.names = c("dlow", "dhigh"))
  # O.data_dhigh_dlow1[]
}

test.model.fits.stratify <- function() {
  require("data.table")
  require("stremr")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)

  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  # OdataNoCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  data(OdataNoCENS)
  OdataNoCENS[OdataNoCENS[,"t"]%in%16,"lastNat1"] <- NA
  # head(OdataNoCENS)
  OdataNoCENS.DT <- as.data.table(OdataNoCENS, key=c("ID", "t"))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS.DT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS.DT[is.na(C), "C" := 0]
  OdataNoCENS.DT[is.na(N), "N" := 0]
  OdataNoCENS.DT[is.na(TI), "TI" := 0]
  OdataNoCENS.DT[is.na(highA1c), "highA1c" := 0]
  OdataNoCENS.DT[is.na(lastNat1), "lastNat1" := 0]
  OdataNoCENS.DT[, "TI_1" := 1L]

  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  res <- stremr(OdataNoCENS.DT, ID = "ID", t_name = "t", tvals = c(0:2),
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          intervened_TRT = "TI_1",
          gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR,
          start_h2o_cluster = FALSE)
  # res$estimates
  # res$dataDT
  # res$wts_data
  # res$OData.R6

  ## ------------------------------------------------------------------------
  ## TEST FOR PROBLEMS WITH > 1 REGRESSION AND > 1 STATA (SHOULD WORK)
  ## ------------------------------------------------------------------------
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  stratify_TRT <- list(
    TI=c("t == 0L", "t == 1L", "t == 2L", "t == 3L", "t == 4L", "t == 5L", "t > 5L")
    )

  gform_CENS_test <- c("C1 ~ highA1c", "C2 ~ highA1c")
  stratify_CENS_test <- list(C1=c("t < 16", "t == 16"), C2=c("t < 16", "t == 16"))
  Odat_DT_test <- OdataNoCENS.DT
  Odat_DT_test[, "C1" := C]
  Odat_DT_test[, "C2" := C]
  OData <- importData(Odat_DT_test, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = c("C1","C2"), TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData, gform_CENS = gform_CENS_test, stratify_CENS = stratify_CENS_test, gform_TRT = gform_TRT,
                                stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)


  # --------------------------------
  # EXAMPLE 2:
  # *********************************************
  # This needs closer look:
  # * We are defining CENS as "C", but then gform_CENS has 3 covars (C, TI, N)
  # * This should either give an error, or we should just drop CENS argument entirely and specify CENS vars based on outcomes of
  # * gform_CENS
  # * Same applies to TRT & MONITOR
  # *********************************************
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify_CENS <- rep(list(strat.str), 3)
  names(stratify_CENS) <- c("C", "TI", "N")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  # options(stremr.verbose = TRUE)
  res <- stremr(OdataNoCENS.DT, ID = "ID", t_name = "t", tvals = c(0:2),
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
                gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
                gform_TRT = gform_TRT,
                gform_MONITOR = gform_MONITOR,
                intervened_TRT = "TI_1",
                # MSMGLMpkg = "h2o",
                tbreaks = c(4,6),
                start_h2o_cluster = FALSE)
                # noCENScat = 0L)

  # res$estimates
  # res$dataDT
  # res$wts_data
  # res$OData.R6

  # --------------------------------
  # EXAMPLE 3:
  # --------------------------------
  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  strat.str <- c("t == 0L", "t > 0")
  stratify_CENS <- rep(list(strat.str), 3)
  names(stratify_CENS) <- c("C", "TI", "N")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  # system.time(
  res <- stremr(OdataNoCENS.DT, ID = "ID", t_name = "t", tvals = c(0:2),
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
                gform_CENS = gform_CENS,
                stratify_CENS = stratify_CENS,
                gform_TRT = gform_TRT,
                gform_MONITOR = gform_MONITOR,
                intervened_TRT = "TI_1",
                tbreaks = c(4,6),
                start_h2o_cluster = FALSE)
  # res$estimates
  # res$dataDT
  # res$wts_data
  # res$OData.R6

  # --------------------------------
  # EXAMPLE 4:
  # --------------------------------
  # define lagged TI, first value is always 1 (always monitored at the first time point):
  OdataNoCENS.DT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]

  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify_CENS <- list(C = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(TI=c("t == 0L", "t > 0L & TI.tminus1 == 0L", "t > 0L & TI.tminus1 == 1L"))
  # **** really want to define it like this ****
  # gform_TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
  #               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))
  gform_MONITOR <- "N ~ 1"

  res <- stremr(OdataNoCENS.DT, ID = "ID", t_name = "t", tvals = c(0:2),
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
                gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
                gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
                gform_MONITOR = gform_MONITOR,
                intervened_TRT = "TI_1",
                tbreaks = c(4,6),
                start_h2o_cluster = FALSE)

  # res$estimates
  # res$dataDT
  # res$wts_data
  # res$OData.R6
}

test.error.fits.stratify <- function() {
  require("data.table")
  require("stremr")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)

  data(OdataNoCENS)
  OdataNoCENS[OdataNoCENS[,"t"]%in%16,"lastNat1"] <- NA
  # head(OdataNoCENS)
  OdataNoCENS.DT <- as.data.table(OdataNoCENS, key=c("ID", "t"))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS.DT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]

  # --------------------------------
  # EXAMPLE 5: Test for error when item names in stratification list do not match the outcome names in regression formula(s)
  # --------------------------------
  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify_CENS <- list(wrongC = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  checkException(
      stremr(OdataNoCENS.DT, ID = "ID", t_name = "t", tvals = c(0:2),
            covars = c("highA1c", "lastNat1"),
            CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
            gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
            gform_TRT = gform_TRT,
            gform_MONITOR = gform_MONITOR,
            intervened_TRT = "TI_1",
            tbreaks = c(4,6),
            start_h2o_cluster = FALSE)
      )

}