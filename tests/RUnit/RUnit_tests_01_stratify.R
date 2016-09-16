simulateDATA.fromDAG <- function(catC = FALSE, Nsize = 1000, rndseed = 124356){
  require("simcausal")
  options(simcausal.verbose = FALSE)
  # get the likelihood under N(t) Bernoulli with P(N(t)=1)=p:
  g.p <- function(p){
    function(O.data, Yname){
      g.N <- rep(p, nrow(O.data))
      g.N <- ifelse(O.data[, "N.t"] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
      return(g.N)
    }
  }
  Drm <- DAG.empty()
  Drm <- Drm +
    node("Y", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  #Z(0)
    node("lastNat1", t = 0, distr = "rconst", const = 0) +  #Z(0) - see below for definition,  set at 0 at t = 0 by convention (see also below)
    node("highA1c.UN", t = 0, distr = "rbern", prob = 0.05) +  # I(0) - possibly unobserved lab
    node("highA1c", t = 0, distr = "rbern", prob = highA1c.UN[0]) +  # I*(0) = I(0)*T(-1) with T(-1) = 1 by convention - all patient come with an A1c value known - T(-1) is what I call N(-1)
    node("CVD", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 1, 0.5, 0.1)) +  #Z(0)
    node("timelowA1c.UN", t = 0, distr = "rbern", prob = 1-highA1c.UN[0]) +  # Z(0) - counts the number of time the (possibly unobserved) A1c was low
    node("TI", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 0, ifelse(CVD[0] == 1, 0.5, 0.1), ifelse(CVD[0] == 1, 0.9, 0.5))) +
    node("CatC", t = 0, distr = "rconst", const = 0) +
    node("C", t = 0, distr = "rconst", const = ifelse(CatC[t] > 0L, 1, 0), EFU = TRUE) +  # last time point is when admin end of study occurs
    node("N", t = 0, distr = "rbern", prob = 1) +
    node("Y", t = 1:16, distr = "rbern", prob = plogis(-6.5 + 1*CVD[0] + 4*highA1c.UN[t-1] + 0.05*timelowA1c.UN[t-1]),  EFU = TRUE) +  # Z(t)
    node("lastNat1", t = 1:16, distr = "rconst", const = ifelse(N[t-1] == 0, lastNat1[t-1] + 1, 0)) +  # Z(1)  just a function of past \bar{N}(t-1) - 0 probs current N at 1,  1 probs previous N a 1,  2 probs  the one before the previous was at 1,  etc.
    node("highA1c.UN", t = 1:16, distr = "rbern", prob = ifelse(TI[t-1] == 1, 0.1, ifelse(highA1c.UN[t-1] == 1, 0.9, min(1, 0.1 + t/16)))) +  # I(t)
    node("highA1c", t = 1:16, distr = "rbern", prob = ifelse(N[t-1] == 1, highA1c.UN[t], highA1c[t-1])) + # I*(m)=I(m)*T(m-1)  (I actually replace I*(m)=0 with when T(m-1)=0 with last value carried forward,  i.e. I*(m-1)
    node("timelowA1c.UN", t=1:16, distr="rnorm", mean=sum(1-highA1c.UN[0:t]),  sd=0) +  # Z(m)
    node("TI", t = 1:16, distr = "rbern",
      prob =
        ifelse(TI[t-1] == 1, 1,
          ifelse(N[t-1] == 1,
            ifelse(highA1c[t] == 0,
              ifelse(CVD[0] == 1, 0.3, 0.1),
                ifelse(CVD[0] == 1, 0.7, 0.5)), 0))) +
    node("CatC", t = 1:16, distr = "rconst", const = 0L) +
    node("C", t = 1:16, distr = "rconst", const = ifelse(CatC[t] > 0L, 1L, 0L), EFU = TRUE) +  # last time point is when admin end of study occurs
###### N(t)
    node("N", t = 1:16, distr = "rbern", prob = 1)
  if (catC) {
    Drm <- Drm + node("CatC", t = 1:16, distr = "rcat.b0", probs = c(0.90, 0.05, 0.05))
  }
  g05.Params <- list(name = "Sporadic 0.5", gInt.N = g.p(0.5), Nprob.t0 = 0.5, Nprob.tplus = 0.5)
  prob.t0 <- g05.Params$Nprob.t0
  prob.tplus <- g05.Params$Nprob.tplus
  if (!is.null(prob.t0)) {
    Drm <- Drm + node("N", t = 0, distr = "rbern", prob = .(prob.t0))
    Drm <- Drm + node("N", t = 1:16, distr = "rbern", prob = .(prob.tplus))
  }
  DAGobj <- set.DAG(Drm, latent.v = c("highA1c.UN", "timelowA1c.UN"))
  O.data <- sim(DAG = DAGobj, n = Nsize, wide = FALSE, rndseed = rndseed)
  O.data <- O.data[,!names(O.data)%in%c("highA1c.UN", "timelowA1c.UN")]
  # head(O.data)

  # -------------------------------------------------------------------------------------------
  # convert lastNat1 to integer
  # -------------------------------------------------------------------------------------------
  O.data <- data.table(O.data)
  O.data[, "lastNat1" := as.integer(lastNat1)]
  # -------------------------------------------------------------------------------------------
  # Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
  # -------------------------------------------------------------------------------------------
  OUTCOME <- "Y"
  shifted.OUTCOME <- "Y.tplus1"
  O.data[, (shifted.OUTCOME) := shift(get(OUTCOME), n = 1L, type = "lead"), by = ID]
  O.data <- O.data[!get(OUTCOME)%in%1,]
  O.data <- O.data[,(OUTCOME) := NULL]
  return(data.frame(O.data))
}

notrun.save.example.data.01 <- function() {
  OdataNoCENS <- simulateDATA.fromDAG(Nsize = 1000, rndseed = 124356)
  head(OdataNoCENS)
  # --------------------------------
  # save data as csv
  # --------------------------------
  # data.table::fwrite(OdataNoCENS, "./OdataNoCENS.csv", turbo = TRUE, verbose = TRUE, na = "NA_h2o")
  # --------------------------------
  # save as compressed R file
  # --------------------------------
  require("tools")
  save(OdataNoCENS, compress = TRUE, file = "./data/OdataNoCENS.rda", compression_level = 9)
  resaveRdaFiles("./data/OdataNoCENS.rda", compress = "bzip2")
}

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
  OdataNoCENS.DT <- as.data.table(OdataNoCENS, key=c(ID, t))
  addN.t1 <- defineMONITORvars(OdataNoCENS, ID = "ID", t = "t", imp.I = "N",
                        MONITOR.name = "N.new", tsinceNis1 = "last.Nt")
  # addN.t1[]
  addN.t2 <- defineMONITORvars(OdataNoCENS.DT, ID = "ID", t = "t", imp.I = "N",
                               MONITOR.name = "N.new", tsinceNis1 = "last.Nt")
  # addN.t2[]
  OdataNoCENS_dhigh_dlow1 <- defineTRTrules(OdataNoCENS, theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
                                      CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1", rule.names = c("dlow", "dhigh"))
  # OdataNoCENS_dhigh_dlow1[]

  OdataNoCENS_dhigh_dlow2 <- defineTRTrules(OdataNoCENS.DT, theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
                                        CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1", rule.names = c("dlow", "dhigh"))
  # OdataNoCENS_dhigh_dlow2[]

  setnames(OdataNoCENS.DT,old = "ID",new = "ID.expression")
  addN.t1 <- defineMONITORvars(OdataNoCENS.DT, ID = "ID.expression", t = "t", imp.I = "N", MONITOR.name = "N.new", tsinceNis1 = "last.Nt")
  # addN.t1[]

  OdataNoCENS_dhigh_dlow1 <- defineTRTrules(OdataNoCENS.DT, theta = c(0,1), ID = "ID.expression", t = "t", I = "highA1c",
                                        CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1", rule.names = c("dlow", "dhigh"))
  # O.data_dhigh_dlow1[]
}

test.model.fits.stratify <- function() {
  require("data.table")
  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  # OdataNoCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  data(OdataNoCENS)
  OdataNoCENS[OdataNoCENS[,"t"]%in%16,"lastNat1"] <- NA
  # head(OdataNoCENS)
  OdataNoCENS.DT <- as.data.table(OdataNoCENS, key=c(ID, t))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS.DT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]

  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  res <- stremr(OdataNoCENS.DT, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
          gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR)
  # res$IPW_estimates
  # res$dataDT
  # res$wts_data
  # res$OData.R6

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
  res <- stremr(OdataNoCENS.DT, ID = "ID", t = "t",
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
                gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
                gform_TRT = gform_TRT,
                gform_MONITOR = gform_MONITOR)
                # noCENScat = 0L)

  # res$IPW_estimates
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
  res <- stremr(OdataNoCENS.DT, ID = "ID", t = "t",
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
                gform_CENS = gform_CENS,
                stratify_CENS = stratify_CENS,
                gform_TRT = gform_TRT,
                gform_MONITOR = gform_MONITOR)
  # res$IPW_estimates
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

  res <- stremr(OdataNoCENS.DT, ID = "ID", t = "t",
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
                gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
                gform_TRT = gform_TRT, stratify_TRT = stratify_TRT,
                gform_MONITOR = gform_MONITOR)

  # res$IPW_estimates
  # res$dataDT
  # res$wts_data
  # res$OData.R6
}

test.error.fits.stratify <- function() {
  require("data.table")
  data(OdataNoCENS)
  OdataNoCENS[OdataNoCENS[,"t"]%in%16,"lastNat1"] <- NA
  # head(OdataNoCENS)
  OdataNoCENS.DT <- as.data.table(OdataNoCENS, key=c(ID, t))
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
      stremr(OdataNoCENS.DT, ID = "ID", t = "t",
            covars = c("highA1c", "lastNat1"),
            CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1",
            gform_CENS = gform_CENS, stratify_CENS = stratify_CENS,
            gform_TRT = gform_TRT,
            gform_MONITOR = gform_MONITOR))

}