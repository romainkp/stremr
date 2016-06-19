library("R6")
library("data.table")

test.buildingblocks <- function() {
  library("data.table")
  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  Nsize <- 1000
  O.data <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  O.data[O.data[,"t"]%in%16,"lastNat1"] <- NA
  O.data <- O.data[,!names(O.data)%in%c("highA1c.UN", "timelowA1c.UN")]
  # head(O.data)
  O.data.DT <- as.data.table(O.data, key=c(ID, t))

  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  OData <- importData(O.data, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
  # have to use this to print:
  # OData$dat.sVar[]
  modelfits.g0 <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  # , gstar_TRT = ..., gstar_MONITOR = ...
  wts.DT <- getIPWeights(modelfits.g0 = modelfits.g0, OData = OData)
  survNP_ests <- get_survNP(wts.DT, OData)
  survNP_ests$IPW_estimates
  # survMSM(data.wts.list, t, MSMregform)

  res <- stremr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
          # noCENScat = 0L)
  res$IPW_estimates
  all.equal(survNP_ests$IPW_estimates, survNP_ests$IPW_estimates)
  # res$dataDT

  # --------------------------------
  # EXAMPLE 2:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify_CENS <- rep(list(strat.str), 3)
  names(stratify_CENS) <- c("C", "TI", "N")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  OData <- importData(O.data, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
  modelfits.g0 <- fitPropensity(OData = OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  wts.DT <- getIPWeights(modelfits.g0 = modelfits.g0, OData = OData)
    # , gstar_TRT = ..., gstar_MONITOR = ...)
  survNP_ests <- get_survNP(wts.DT, OData)
  # survMSM(data.wts.list, t, MSMregform)
  survNP_ests$IPW_estimates

  # system.time(
  # res <-
  #   stremr(O.data, ID = "ID", t = "t",
  #         covars = c("highA1c", "lastNat1"),
  #         CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
  #         gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
  #         gform_TRT = gform_TRT,
  #         gform_MONITOR = gform_MONITOR)
  #         # noCENScat = 0L)
  #   )
  # res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 3:
  # --------------------------------
  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify_CENS <- rep(list(strat.str), 3)
  names(stratify_CENS) <- c("C", "TI", "N")
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"

  OData <- importData(O.data, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
  modelfits.g0 <- fitPropensity(OData = OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  wts.DT <- getIPWeights(modelfits.g0 = modelfits.g0, OData = OData)
    # , gstar_TRT = ..., gstar_MONITOR = ...)
  survNP_ests <- get_survNP(wts.DT, OData)
  # survMSM(data.wts.list, t, MSMregform)
  survNP_ests$IPW_estimates

  # system.time(
  # res <-
  #   stremr(O.data, ID = "ID", t = "t",
  #         covars = c("highA1c", "lastNat1"),
  #         CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
  #         gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
  #         gform_TRT = gform_TRT,
  #         gform_MONITOR = gform_MONITOR)
  #         # noCENScat = 0L)
  #   )
  # res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 4:
  # --------------------------------
  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify_CENS <- list(C = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(TI=c("t == 0L", "t > 0L & TI.tminus1 == 0L", "t > 0L & TI.tminus1 == 1L"))
  # **** really want to define it like this ****
  # gform_TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
  #               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))
  gform_MONITOR <- "N ~ 1"

  system.time(
  res <-
    stremr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
          gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
          gform_MONITOR = gform_MONITOR)
          # noCENScat = 0L)
    )
  res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 5: Test for error when item names in stratification list do not match the outcome names in regression formula(s)
  # --------------------------------
  gform_CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify_CENS <- list(wrongC = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  checkException(
      stremr(O.data, ID = "ID", t = "t",
            covars = c("highA1c", "lastNat1"),
            CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
            gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
            gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
            gform_MONITOR = gform_MONITOR))
}
