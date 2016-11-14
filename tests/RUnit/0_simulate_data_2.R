
# ------------------------------------------------------------------------------------------------------
# SIMULATE, SAVE AND COMPRESS THE DATASET FROM THE EXAMPLE
# ------------------------------------------------------------------------------------------------------
simulate.sim.DATA.fromDAG.prepDATA <- function(Nsize = 10000, rndseed = 55466) {
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  require("simcausal")
  # -----------------------------------------------------------
  # SIMULATION PARAMS:
  # -----------------------------------------------------------
  prob.t0 <- prob.tplus <- 0.5
  ###########################################################################################################
  # Run the simulation for one scenario on N(t) and perform estimation for two scenarios:
  ###########################################################################################################
  set.DAG.rm.template <- function(){
    require("simcausal")
    options(simcausal.verbose = FALSE)
    Drm <- DAG.empty()
    Drm <- Drm +
      node("Y", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  #Z(0)
      node("lastNat1", t = 0, distr = "rconst", const = 0) +  #Z(0) - see below for definition,  set at 0 at t = 0 by convention (see also below)
      node("highA1c.UN", t = 0, distr = "rbern", prob = 0.05) +  # I(0) - possibly unobserved lab
      node("highA1c", t = 0, distr = "rbern", prob = highA1c.UN[0]) +  # I*(0) = I(0)*T(-1) with T(-1) = 1 by convention - all patient come with an A1c value known - T(-1) is what I call N(-1)
      node("CVD", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 1, 0.5, 0.1)) +  #Z(0)
      node("timelowA1c.UN", t = 0, distr = "rbern", prob = 1-highA1c.UN[0]) +  # Z(0) - counts the number of time the (possibly unobserved) A1c was low
      node("TI", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 0, ifelse(CVD[0] == 1, 0.5, 0.1), ifelse(CVD[0] == 1, 0.9, 0.5))) +
      node("C", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  # no censoring in first bin
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
      node("C", t = 1:16, distr = "rbern", prob = ifelse(t == 16, 1, 0), EFU = TRUE) +  # last time point is when admin end of study occurs
      node("N", t = 1:16, distr = "rbern", prob = 1)
    return(Drm)
  }

  DAGrm <- set.DAG.rm.template()
  vecfun.add("N.Pois.1"); vecfun.add("N.Pois.3"); vecfun.add("N.sporadic.2"); vecfun.add("N.sporadic.4")
  #### Simulate with P(N(0)=1) = prob.t0 and P(N(t>0)=1) = prob.tplus
  if (!is.null(prob.t0)) {
    DAGrm <- DAGrm + node("N", t = 0, distr = "rbern", prob = .(prob.t0))
    DAGrm <- DAGrm + node("N", t = 1:16, distr = "rbern", prob = .(prob.tplus))
  }
  DAGrm <- set.DAG(DAGrm)
  Odat <- sim(DAG = DAGrm, n = Nsize, wide = FALSE, rndseed = rndseed)
  Odat[Odat[,"t"]%in%16,"lastNat1"] <- NA
  Odat <- Odat[,!names(Odat)%in%c("highA1c.UN", "timelowA1c.UN")]

  # -------------------------------------------------------------------------------------------
  # make a data.table
  # -------------------------------------------------------------------------------------------
  ID <- "ID"; t <- "t"; TRT <- "TI"; CENS <- "C"; MONITOR <- "N"; outcome <- "Y"; I <- "highA1c";
  OdatDT <- data.table(Odat, key = c(ID, t))

  OdatDT[, "lastNat1" := as.integer(lastNat1)] # convert to integer
  # -------------------------------------------------------------------------------------------
  # Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
  # -------------------------------------------------------------------------------------------
  OUTCOME <- "Y"
  shifted.OUTCOME <- OUTCOME%+%".tplus1"
  OdatDT[, (shifted.OUTCOME) := shift(get(OUTCOME), n = 1L, type = "lead"), by = ID]
  OdatDT <- OdatDT[!get(OUTCOME)%in%1,]
  OdatDT <- OdatDT[,(OUTCOME) := NULL]
  # setnames(OdatDT,old = shifted.OUTCOME, new = OUTCOME)

  # --------------------------------
  # Define two dynamic regimes (counterfactual treatment assignment under two rules (dlow & dhigh))
  # --------------------------------
  # Counterfactual TRT assignment for rule dlow (equivalent to always treated):
  rule_name1 <- "dlow"
  OdatDT[,"gTI." %+% rule_name1 := 1L]
  # Counterfactual TRT assignment for dynamic rule dhigh -> start TRT only when I=1 (highA1c = 1)
  rule_name2 <- "dhigh"
  OdatDT_TIdhigh <- stremr::defineIntervedTRT(OdatDT, theta = 1, ID = ID,
                                            t = t, I = I, CENS = CENS, TRT = TRT,
                                            MONITOR = MONITOR,
                                            tsinceNis1 = "lastNat1",
                                            new.TRT.names = "gTI." %+% rule_name2)
  OdatDT <- merge(OdatDT, OdatDT_TIdhigh, by=c(ID, t))

  # --------------------------------
  # Define monitoring regimen probabilit(ies):
  # --------------------------------
  # N^*(t) Bernoulli with P(N^*(t)=1)=p
  g.p <- function(Odat, p) return(rep(p, nrow(Odat)))
  # N^*(t) Poisson:
  g.Pois <- function(Odat, lambda, lastNat1 = "lastNat1") {
    g.N <- 1 / (ppois(Odat[[lastNat1]], lambda, lower.tail = FALSE) / dpois(Odat[[lastNat1]], lambda) + 1)
    g.N[is.na(g.N)] <- 0
    return(g.N)
  }
  Odat_DT <- OdatDT[, c("gPois3.yrly", "gPois3.biyrly", "gp05") := list(g.Pois(OdatDT, lambda = 3), g.Pois(OdatDT, lambda = 1), g.p(OdatDT, p = 0.5))][]

  return(Odat_DT)
}

# ------------------------------------------------------------------------------------------------------------------
# run this to resave the simulated dataset as a data.table (with all post-processing variables defined for stremr)
# ------------------------------------------------------------------------------------------------------------------
notrun.save.example.data.04 <- function() {
  OdatDT_10K <- simulate.sim.DATA.fromDAG.prepDATA(Nsize = 10000, rndseed = 55466)
  # --------------------------------
  # save data as csv
  # --------------------------------
  # data.table::fwrite(OdatDT_10K, "./OdatDT_10K.csv", verbose = TRUE, na = "NA_h2o")
  # --------------------------------
  # save as compressed R file
  # --------------------------------
  require("tools")
  save(OdatDT_10K, compress = TRUE, file = "./data/OdatDT_10K.rda", compression_level = 9)
  resaveRdaFiles("./data/OdatDT_10K.rda", compress = "bzip2")
}
