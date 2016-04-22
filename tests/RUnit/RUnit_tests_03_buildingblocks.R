library("R6")
library("data.table")

# DatNetTEMP <- R6Class(classname = "DatNetTEMP",
#   portable = TRUE,
#   class = TRUE,
#   public = list(
#     dataDT = NULL,
#     initialize = function() {
#       invisible(self)
#     }
#   ),
#   active = list(
#     temp = function() { NULL }
#   ),
#   private = list(
#     .tmp = NULL
#   )
# )

# returnDT <- function() {
#   DT <- data.table(A = c(1,2,3), B = c("a", "b", "c"))
#   DatNetTEMP$dataDT <- DT
#   return(list(dataDT = DT, DatNet = DatNetTEMP))
#   # return(list(DT = DT, rn = runif(100)))
# }

# datOBJR6 <- returnDT()
# datOBJR6$dataDT
# datOBJR6$DatNet$dataDT
# datDT$DT

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
  gform.CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"
  OData <- get_Odata(O.data, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
  # have to use this to print:
  # OData$dat.sVar[]
  modelfits.g0 <- get_fits(OData = OData, gform.CENS = gform.CENS, gform.TRT = gform.TRT, gform.MONITOR = gform.MONITOR)
  # , gstar.TRT = ..., gstar.MONITOR = ...
  wts.DT <- get_weights(modelfits.g0 = modelfits.g0, OData = OData)
  survNP_ests <- get_survNP(wts.DT, OData)
  survNP_ests$IPW_estimates
  # get_survMSM(data.wts.list, t, MSMregform)

  res <-
    stremr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform.CENS = gform.CENS, gform.TRT = gform.TRT, gform.MONITOR = gform.MONITOR)
          # noCENS.cat = 0L)
  res$IPW_estimates
  all.equal(survNP_ests$IPW_estimates, survNP_ests$IPW_estimates)
  # res$dataDT

  # --------------------------------
  # EXAMPLE 2:
  # --------------------------------
  gform.CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify.CENS <- rep(list(strat.str), 3)
  names(stratify.CENS) <- c("C", "TI", "N")
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"

  OData <- get_Odata(O.data, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
  modelfits.g0 <- get_fits(OData = OData, gform.CENS = gform.CENS, stratify.CENS = stratify.CENS, gform.TRT = gform.TRT, gform.MONITOR = gform.MONITOR)
  get_weights(modelfits.g0 = modelfits.g0, OData = OData, gstar.TRT = ..., gstar.MONITOR = ...)
  get_survNP(data.wts, t)
  get_survMSM(data.wts.list, t, MSMregform)

  # system.time(
  # res <-
  #   stremr(O.data, ID = "ID", t = "t",
  #         covars = c("highA1c", "lastNat1"),
  #         CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
  #         gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
  #         gform.TRT = gform.TRT,
  #         gform.MONITOR = gform.MONITOR)
  #         # noCENS.cat = 0L)
  #   )
  # res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 3:
  # --------------------------------
  gform.CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"

  gform.CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify.CENS <- rep(list(strat.str), 3)
  names(stratify.CENS) <- c("C", "TI", "N")
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"

  OData <- get_Odata(O.data, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
  modelfits.g0 <- get_fits(OData = OData, gform.CENS = gform.CENS, stratify.CENS = stratify.CENS, gform.TRT = gform.TRT, gform.MONITOR = gform.MONITOR)
  get_weights(modelfits.g0 = modelfits.g0, OData = OData, gstar.TRT = ..., gstar.MONITOR = ...)
  get_survNP(data.wts, t)
  get_survMSM(data.wts.list, t, MSMregform)

  # system.time(
  # res <-
  #   stremr(O.data, ID = "ID", t = "t",
  #         covars = c("highA1c", "lastNat1"),
  #         CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
  #         gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
  #         gform.TRT = gform.TRT,
  #         gform.MONITOR = gform.MONITOR)
  #         # noCENS.cat = 0L)
  #   )
  # res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 4:
  # --------------------------------
  gform.CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify.CENS <- list(C = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  stratify.TRT <- list(TI=c("t == 0L", "t > 0L & TI.tminus1 == 0L", "t > 0L & TI.tminus1 == 1L"))
  # **** really want to define it like this ****
  # gform.TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
  #               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))
  gform.MONITOR <- "N ~ 1"

  system.time(
  res <-
    stremr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
          gform.TRT = gform.TRT, stratify.TRT = stratify.TRT,
          gform.MONITOR = gform.MONITOR)
          # noCENS.cat = 0L)
    )
  res$IPW_estimates
  # res$dataDT

  # --------------------------------
  # EXAMPLE 5: Test for error when item names in stratification list do not match the outcome names in regression formula(s)
  # --------------------------------
  gform.CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  stratify.CENS <- list(wrongC = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
  checkException(
      stremr(O.data, ID = "ID", t = "t",
            covars = c("highA1c", "lastNat1"),
            CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
            gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
            gform.TRT = gform.TRT, stratify.TRT = stratify.TRT,
            gform.MONITOR = gform.MONITOR))
}
