set.DAG.rm.template <- function(catC = FALSE){
  require(simcausal)
  options(simcausal.verbose = FALSE)
  rcategor.int.b0 <- function(n, probs) { rcategor.int(n, probs)-1 }
  Drm <- DAG.empty()
###### L(0)=(Z(0),I*(0)) - I*(m) represents value of test at m-1 if test performed: I*(m)=I(m)*T(m-1) where I(m) is value of test (possibly unobserved)
  Drm <- Drm +
    node("Y", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  #Z(0)
    node("lastNat1", t = 0, distr = "rconst", const = 0) +  #Z(0) - see below for definition,  set at 0 at t = 0 by convention (see also below)
    node("highA1c.UN", t = 0, distr = "rbern", prob = 0.05) +  # I(0) - possibly unobserved lab
    node("highA1c", t = 0, distr = "rbern", prob = highA1c.UN[0]) +  # I*(0) = I(0)*T(-1) with T(-1) = 1 by convention - all patient come with an A1c value known - T(-1) is what I call N(-1)
    node("CVD", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 1, 0.5, 0.1)) +  #Z(0)
    node("timelowA1c.UN", t = 0, distr = "rbern", prob = 1-highA1c.UN[0]) +  # Z(0) - counts the number of time the (possibly unobserved) A1c was low
###### A(0)
    node("TI", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 0, ifelse(CVD[0] == 1, 0.5, 0.1), ifelse(CVD[0] == 1, 0.9, 0.5))) +
    node("C", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  # no censoring in first bin
###### T(0)
    node("N", t = 0, distr = "rbern", prob = 1) +
###### L(t) for t = 1, ..., 16
    node("Y", t = 1:16, distr = "rbern", prob = plogis(-6.5 + 1*CVD[0] + 4*highA1c.UN[t-1] + 0.05*timelowA1c.UN[t-1]),  EFU = TRUE) +  # Z(t)
    node("lastNat1", t = 1:16, distr = "rconst", const = ifelse(N[t-1] == 0, lastNat1[t-1] + 1, 0)) +  # Z(1)  just a function of past \bar{N}(t-1) - 0 probs current N at 1,  1 probs previous N a 1,  2 probs  the one before the previous was at 1,  etc.
    node("highA1c.UN", t = 1:16, distr = "rbern", prob = ifelse(TI[t-1] == 1, 0.1, ifelse(highA1c.UN[t-1] == 1, 0.9, min(1, 0.1 + t/16)))) +  # I(t)
    node("highA1c", t = 1:16, distr = "rbern", prob = ifelse(N[t-1] == 1, highA1c.UN[t], highA1c[t-1])) + # I*(m)=I(m)*T(m-1)  (I actually replace I*(m)=0 with when T(m-1)=0 with last value carried forward,  i.e. I*(m-1)
    node("timelowA1c.UN", t=1:16, distr="rnorm", mean=sum(1-highA1c.UN[0:t]),  sd=0) +  # Z(m)
###### A(t)
    node("TI", t = 1:16, distr = "rbern",
      prob =
        ifelse(TI[t-1] == 1, 1,
          ifelse(N[t-1] == 1,
            ifelse(highA1c[t] == 0,
              ifelse(CVD[0] == 1, 0.3, 0.1),
                ifelse(CVD[0] == 1, 0.7, 0.5)), 0))) +
    node("C", t = 1:16, distr = "rbern", prob = ifelse(t == 16, 1, 0), EFU = TRUE) +  # last time point is when admin end of study occurs
###### N(t)
    node("N", t = 1:16, distr = "rbern", prob = 1)

  if (catC) {
    Drm <- Drm + node("C", t = 1:16, distr = "rcategor.int.b0", probs = c(0.95, 0.025, 0.025), EFU = TRUE)
  }
  return(Drm)
}

test.model.fits.categorCENSOR <- function() {
  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  Nsize <- 50000
  library("simcausal")
  DAGrm <- set.DAG.rm.template(catC=TRUE)

  g05.Params = list(name = "Sporadic 0.5", gInt.N = g.p(0.5), Nprob.t0 = 0.5, Nprob.tplus = 0.5)
  prob.t0 <- g05.Params$Nprob.t0
  prob.tplus <- g05.Params$Nprob.tplus
  if (!is.null(prob.t0)) {
    DAGrm <- DAGrm + node("N", t = 0, distr = "rbern", prob = .(prob.t0))
    DAGrm <- DAGrm + node("N", t = 1:16, distr = "rbern", prob = .(prob.tplus))
  }
  DAG <- set.DAG(DAGrm)
  O.data <- sim(DAG = DAG, n = Nsize, wide = FALSE, rndseed = 124356)
  O.data[O.data[,"t"]%in%16,"lastNat1"] <- NA
  O.data <- O.data[,!names(O.data)%in%c("highA1c.UN", "timelowA1c.UN")]
  head(O.data)
  unique(O.data[,"C"])
  # [1]  0  2  1 NA

  # --------------------------------
  # EXAMPLE 1:
  # --------------------------------
  gform.CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"
  system.time(
  res <-
    estimtr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform.CENS = gform.CENS, gform.TRT = gform.TRT, gform.MONITOR = gform.MONITOR)
          # noCENS.cat = 0L)
    )
  summary(res$h_gN)
  #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
  # 0.0008805 0.0405500 0.1581000 0.1575000 0.2271000 1.0000000

  # --------------------------------
  # EXAMPLE 2:
  # --------------------------------
  gform.CENS <- "C + TI + N ~ highA1c + lastNat1"
  strat.str <- c("t == 0L", "t > 0")
  stratify.CENS <- rep(list(strat.str), 3)
  names(stratify.CENS) <- c("C", "TI", "N")
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"

  system.time(
  res <-
    estimtr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
          gform.TRT = gform.TRT,
          gform.MONITOR = gform.MONITOR)
          # noCENS.cat = 0L)
    )
  summary(res$h_gN)
  #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
  # 0.001018 0.038620 0.164500 0.168500 0.254200 1.000000

  # --------------------------------
  # EXAMPLE 3:
  # --------------------------------
  gform.CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
  gform.TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform.MONITOR <- "N ~ 1"
  system.time(
  res <-
    estimtr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
          gform.TRT = gform.TRT,
          gform.MONITOR = gform.MONITOR)
          # noCENS.cat = 0L)
    )
  summary(res$h_gN)
  #     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
  # 0.001018 0.038620 0.164500 0.168500 0.254200 1.000000

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
    estimtr(O.data, ID = "ID", t = "t",
          covars = c("highA1c", "lastNat1"),
          CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
          gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
          gform.TRT = gform.TRT, stratify.TRT = stratify.TRT,
          gform.MONITOR = gform.MONITOR)
          # noCENS.cat = 0L)
    )
  summary(res$h_gN)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00109 0.08220 0.28870 0.23970 0.35340 1.00000

}

















