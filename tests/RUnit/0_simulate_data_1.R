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
    node("Y",             t = 0,    distr = "rbern",  prob  = 0, EFU = TRUE) +  #Z(0)
    node("lastNat1",      t = 0,    distr = "rconst", const = 0) +  #Z(0) - see below for definition,  set at 0 at t = 0 by convention (see also below)
    node("highA1c.UN",    t = 0,    distr = "rbern",  prob  = 0.05) +  # I(0) - possibly unobserved lab
    node("highA1c",       t = 0,    distr = "rbern",  prob  = highA1c.UN[0]) +  # I*(0) = I(0)*T(-1) with T(-1) = 1 by convention - all patient come with an A1c value known - T(-1) is what I call N(-1)
    node("CVD",           t = 0,    distr = "rbern",  prob  = ifelse(highA1c[0] == 1, 0.5, 0.1)) +  #Z(0)
    node("timelowA1c.UN", t = 0,    distr = "rbern",  prob  = 1-highA1c.UN[0]) +  # Z(0) - counts the number of time the (possibly unobserved) A1c was low
    node("TI",            t = 0,    distr = "rbern",  prob  = ifelse(highA1c[0] == 0, ifelse(CVD[0] == 1, 0.5, 0.1), ifelse(CVD[0] == 1, 0.9, 0.5))) +
    node("CatC",          t = 0,    distr = "rconst", const = 0) +
    node("C",             t = 0,    distr = "rconst", const = ifelse(CatC[t] > 0L, 1, 0), EFU = TRUE) +  # last time point is when admin end of study occurs
    node("N",             t = 0,    distr = "rbern",  prob  = 1) +
    node("Y",             t = 1:16, distr = "rbern",  prob  = plogis(-6.5 + 1*CVD[0] + 4*highA1c.UN[t-1] + 0.05*timelowA1c.UN[t-1]),  EFU = TRUE) +  # Z(t)
    node("lastNat1",      t = 1:16, distr = "rconst", const = ifelse(N[t-1] == 0, lastNat1[t-1] + 1, 0)) +  # Z(1)  just a function of past \bar{N}(t-1) - 0 probs current N at 1,  1 probs previous N a 1,  2 probs  the one before the previous was at 1,  etc.
    node("highA1c.UN",    t = 1:16, distr = "rbern",  prob  = ifelse(TI[t-1] == 1, 0.1, ifelse(highA1c.UN[t-1] == 1, 0.9, min(1, 0.1 + t/16)))) +  # I(t)
    node("highA1c",       t = 1:16, distr = "rbern",  prob  = ifelse(N[t-1] == 1, highA1c.UN[t], highA1c[t-1])) + # I*(m)=I(m)*T(m-1)  (I actually replace I*(m)=0 with when T(m-1)=0 with last value carried forward,  i.e. I*(m-1)
    node("timelowA1c.UN", t=1:16,   distr = "rnorm",  mean  = sum(1-highA1c.UN[0:t]),  sd=0) +  # Z(m)
    node("TI",            t = 1:16, distr = "rbern",
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
  OdataNoCENS <- OdataNoCENS[, !names(OdataNoCENS) %in% "CatC"]
  # --------------------------------
  # save data as csv
  # --------------------------------
  # data.table::fwrite(OdataNoCENS, "./OdataNoCENS.csv", verbose = TRUE, na = "NA_h2o")
  # --------------------------------
  # save as compressed R file
  # --------------------------------
  require("tools")
  save(OdataNoCENS, compress = TRUE, file = "./data/OdataNoCENS.rda", compression_level = 9)
  resaveRdaFiles("./data/OdataNoCENS.rda", compress = "bzip2")
}

# ------------------------------------------------------------------------------------------------------
# SIMULATE, SAVE AND COMPRESS THE DATASET FROM THE EXAMPLE
# ------------------------------------------------------------------------------------------------------
notrun.save.example.data.02 <- function() {
  require("simcausal")
  # Nsize <- 50000
  Nsize <- 1000
  OdataCatCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356, catC=TRUE)
  OdataCatCENS[OdataCatCENS[,"t"]%in%16,"lastNat1"] <- NA
  save(OdataCatCENS, compress = TRUE, file = "./data/OdataCatCENS.rda", compression_level = 9)
  require("tools")
  resaveRdaFiles("./data/OdataCatCENS.rda", compress = "bzip2")

  # --------------------------------
  # save data as csv
  # --------------------------------
  # data.table::fwrite(OdataCatCENS, "./OdataCatCENS.csv", verbose = TRUE, na = "NA_h2o")
}
