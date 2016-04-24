require("data.table")

delta.shift.minus <- seq(-0.9, -0.1, by=0.1)
delta.shift.plus <- seq(0.1, 0.9, by=0.1)
delta.shifts <- c(delta.shift.minus, delta.shift.plus)
# delta.shifts <- c(-0.5,0.5)

# get the likelihood under N(t) Poisson:
g.Pois <- function(O.dataDT, lambda, ID = "ID", OUTCOME = "Y", lastNat1 = "lastNat1", MONITOR = "N"){
  g.N <- 1 / (ppois(O.dataDT[[lastNat1]], lambda, lower.tail = FALSE) / dpois(O.dataDT[[lastNat1]], lambda) + 1)
  g.N <- ifelse(O.dataDT[[MONITOR]] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
  return(g.N)
}
# # get the likelihood under N(t) Poisson:
# g.Pois <- function(lambda){
#   function(O.data, Yname){
#     g.N <- 1 / (ppois(O.data[, "lastN.t"], lambda, lower.tail = FALSE) / dpois(O.data[, "lastN.t"], lambda) + 1)
#     g.N <- ifelse(O.data[,"N.t"] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
#     return(g.N)
#   }
# }
# get the likelihood under N(t) Bernoulli with P(N(t)=1)=p:
g.p <- function(O.dataDT, p, ID = "ID", OUTCOME = "Y", lastNat1 = "lastNat1", MONITOR = "N"){
  g.N <- rep(p, nrow(O.dataDT))
  g.N <- ifelse(O.dataDT[[MONITOR]] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
  return(g.N)
}
# get the likelihood under N(t) Bernoulli with P(N(t)=1)=p:
# g.p <- function(p){
#   function(O.data, Yname){
#     g.N <- rep(p, nrow(O.data))
#     g.N <- ifelse(O.data[, "N.t"] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
#     return(g.N)
#   }
# }
# DAG equation for defining the distribution of N(t) as Poisson (x is time since last visit)
N.Pois.3 <- function(x,lambda = 3){
  g.N <- 1 / (ppois(x, lambda, lower.tail = FALSE) / dpois(x, lambda[1])+1) # will be NA if x=NA
  return(g.N)
}
# DAG equation for defining the distribution of N(t) as Poisson (x is time since last visit)
N.Pois.1 <- function(x,lambda=1){
  g.N <- 1 / (ppois(x, lambda, lower.tail = FALSE) / dpois(x, lambda[1]) + 1) # will be NA if x=NA
  return(g.N)
}

Sim.and.Est.params <- list(
  gcont       = list(name = "Continuous",         gInt.N = g.p(1),        Nprob.t0 = 1,                                  Nprob.tplus = 1,                                     Nt.form = "N.t ~ 1",                  delta.shift = NULL),
  g05         = list(name = "Sporadic 0.5",       gInt.N = g.p(0.5),      Nprob.t0 = 0.5,                                Nprob.tplus = 0.5,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  g035        = list(name = "Sporadic 0.35",      gInt.N = g.p(0.35),     Nprob.t0 = 0.35,                               Nprob.tplus = 0.35,                                  Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  g02         = list(name = "Sporadic 0.2",       gInt.N = g.p(0.2),      Nprob.t0 = 0.2,                                Nprob.tplus = 0.2,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  g01         = list(name = "Sporadic 0.1",       gInt.N = g.p(0.1),      Nprob.t0 = 0.1,                                Nprob.tplus = 0.1,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  # gSpoYear    = list(name = "Sporadic yearly",    gInt.N = g.sporadic.4,  Nprob.t0 = 0,                                  Nprob.tplus = substitute(N.sporadic.4(lastNat1[t])), Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = delta.shifts),
  # gSpoBiyear  = list(name = "Sporadic bi-yearly", gInt.N = g.sporadic.2,  Nprob.t0 = 0.25,                               Nprob.tplus = substitute(N.sporadic.2(lastNat1[t])), Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = 0.1),
  gPois3      = list(name = "Poisson yearly",     gInt.N = g.Pois(3),     Nprob.t0 = substitute(N.Pois.3(lastNat1[t])),  Nprob.tplus = substitute(N.Pois.3(lastNat1[t])),     Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = delta.shifts, gradual.t.vec = 4),
  gPois1      = list(name = "Poisson bi-yearly",  gInt.N = g.Pois(1),     Nprob.t0 = substitute(N.Pois.1(lastNat1[t])),  Nprob.tplus = substitute(N.Pois.1(lastNat1[t])),     Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = 0.1)
  # Two interventions below are never used for generating obs. data, only as g^* for estimating \psi_0:
  # g0delta.shift.p = list(name = "Delta shift p",        gInt.N = g.delta.p,               Nprob.t0 = NULL, Nprob.tplus = NULL, Nt.form = "not used", NOT.USE.AS.g0 = TRUE),
  # g0gradual.t     = list(name = "Shift gradual t",      gInt.N = g.delta.p.byt,           Nprob.t0 = NULL, Nprob.tplus = NULL, Nt.form = "not used", NOT.USE.AS.g0 = TRUE),
  # g0forcevis.t.1  = list(name = "gN.0 visit every t.1", gInt.N = g.delta.forcevisit.t(1), Nprob.t0 = NULL, Nprob.tplus = NULL, Nt.form = "not used", NOT.USE.AS.g0 = TRUE),
  # g0forcevis.t.4  = list(name = "gN.0 visit every t.4", gInt.N = g.delta.forcevisit.t(4), Nprob.t0 = NULL, Nprob.tplus = NULL, Nt.form = "not used", NOT.USE.AS.g0 = TRUE)
)


###########################################################################################################
# Run the simulation for one scenario on N(t) and perform estimation for two scenarios:
###########################################################################################################
# options(width=140)
# to install from github:
# devtools::install_github('osofr/simcausal', build_vignettes = FALSE)
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# require(stremr)
require("simcausal")
require("data.table")
`%+%` <- function(a, b) paste0(a, b)

# wdir <- "/Users/olegsofrygin/Dropbox/KP/monitoring_simstudy/stremr_legacyRcode"
# setwd(wdir)
# source("666_fit_g.C.A.N_calcIPAW_simfcts.R") # estimation functions + DAG for the data generating distribution
# source("666_g_library.R") # scenarios for N(t)
# datadir <- "./"

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

# ---------------------------------------------------------------------------------------------------------
# Get indicators of rule followers for two rules (dlow, dhigh) based on binary I:
# ---------------------------------------------------------------------------------------------------------
add.gtheta.to.O <- function(O.data, ID = "ID", t = "t", TRT = "TI", CENS = "C", MONITOR = "N", I = "highA1c", rule_name1 = "dlow", rule_name2 = "dhigh"){
  ## Count how many patients follow each rule; follow rule dlow - start TI at low A1c threshold, i.e. start at study entry.
  O.dataDT <- data.table(O.data, key = c(ID, t))
  # Add indicator for rule (dlow) -> start TI when highA1c>=0 and stay on treatment ALL the time (i.e., start at study entry)
  # equivalent to counterfactual regiment \barA=1
  rule_name1 <- "dlow"
  O.dataDT[,(rule_name1) := as.integer(all(get(TRT) %in% c(1,NA))), by=eval(ID)]
  # set last row to 0, since its always either Y=1 or C=1
  O.dataDT[O.dataDT[,.I[.N], by = eval(ID)][["V1"]],(rule_name1) := 0L]
  # set All NAs for treatment (TRT) to be NA for the rule as well
  O.dataDT[O.dataDT[,.I[which(is.na(get(TRT)))], by = eval(ID)][["V1"]],(rule_name1) := NA]
  # Add indicator for rule (dhigh) -> start TRT only when I=1 (highA1c)
  rule_name2 <- "dhigh"
  O.dataDT_dhigh <- stremr::follow.rule.d.DT(O.dataDT, theta = 1, ID = ID,
                                              t = t, I = I, CENS = CENS, TRT = TRT,
                                              MONITOR = MONITOR, rule.names = rule_name2)
  O.dataDT_dhigh[, (rule_name2) := as.integer(get(rule_name2))]
  O.dataDT <- merge(O.dataDT, O.dataDT_dhigh, by=c(ID, t))
  return(O.dataDT)
}
# --------------------------------
# (II) Define event indicator and right-censoring indicator, define total follow-up length
# --------------------------------
define_indicators <- function(O.data, ID = "ID", t = "t", TRT = "TI", CENS = "C", MONITOR = "N", I = "highA1c") {
  # Add indicator Delta=I(T<=C) (Y is 1 at some point for the subject):
  O.dataDT <- data.table(O.data, key = c(ID, t))
  EVENT_IND <- "Delta";
  O.dataDT[,(EVENT_IND) := as.integer(any(get(outcome) %in% 1)), by=eval(ID)]
  # Add indicator Delta=I(T>C) (subject was right-censored at some point):
  CENS_IND <- "AnyCensored"; noCENS.cat <- 0L; CENS <- c("C")
  O.dataDT[, (CENS_IND) := FALSE, by=eval(ID)]
  for (Cvar in CENS) {
    O.dataDT[, (CENS_IND) := get(CENS_IND) | any(!get(Cvar) %in% c(eval(noCENS.cat),NA)), by = eval(ID)]
  }
  O.dataDT[, (CENS_IND) := as.integer(get(CENS_IND))]
  ## Add variable that indicates if TI initiation occured previously, relevant for model of TI continuation
  TIcovarname <- "barTIm1eq0"
  O.dataDT[, (TIcovarname) := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  # O.dataDT[1:100,]
  return(O.dataDT)
}

# -----------------------------------------------------------
# SIMULATION PARAMS:
# -----------------------------------------------------------
Nsize <- 50000
# set to TRUE to only run scenarios dealing with bias for g.N w and w/out past N(t), no intervention on N(t)
run_only_bias_noN <- FALSE
NoIntNt <- FALSE
SimParams <- Sim.and.Est.params
# override the scenario name to run for just one scenario:
scen.names <- "g05"
# scen.names <- c("g01_gcont", "g01_gcont_rndHighA1c", "g01_gcont_rndHighA1cNULL", "g01_gcont_Nstar")
# scen.names <- "gSpoBiyear"
# scen.names <- c("gPois3", "gPois1")
# scen.names <- c("g05", "gPois3")
# scen.names <- c("gSpoYear", "gPois3")

# ------------------------------------------------------------------------------------------------------
# (I) Simulate data with binary I (highA1c)
# ------------------------------------------------------------------------------------------------------
prob.t0 <- SimParams[[scen.names[1]]]$Nprob.t0
prob.tplus <- SimParams[[scen.names[1]]]$Nprob.tplus
print("scen.name: " %+% scen.names[1]); print("prob.t0"); print(prob.t0); print("prob.tplus"); print(prob.tplus)
DAGrm <- set.DAG.rm.template()
vecfun.add("N.Pois.1"); vecfun.add("N.Pois.3"); vecfun.add("N.sporadic.2"); vecfun.add("N.sporadic.4")

#### Simulate with P(N(0)=1) = prob.t0 and P(N(t>0)=1) = prob.tplus
if (!is.null(prob.t0)) {
  DAGrm <- DAGrm + node("N", t = 0, distr = "rbern", prob = .(prob.t0))
  DAGrm <- DAGrm + node("N", t = 1:16, distr = "rbern", prob = .(prob.tplus))
}
DAGrm <- set.DAG(DAGrm)
O.data.simstudy <- sim(DAG = DAGrm, n = Nsize, wide = FALSE, rndseed = 55466)
O.data.simstudy[O.data.simstudy[,"t"]%in%16,"lastNat1"] <- NA
O.data <- O.data.simstudy[,!names(O.data.simstudy)%in%c("highA1c.UN", "timelowA1c.UN")]

ID <- "ID"; t <- "t"; TRT <- "TI"; outcome <- "Y"; I <- "highA1c";

# --------------------------------
# (II) Define event indicator and right-censoring indicator, define total follow-up length
# --------------------------------
O.dataDT <- define_indicators(O.data)
# --------------------------------
# (III) Define rules for dlow & dhigh
# --------------------------------
O.dataDT_rules <- add.gtheta.to.O(O.dataDT)
# --------------------------------
# (IV) Define monitoring regimen probabilit(ies):
# --------------------------------
gstar1.N.Pois3.yearly <- g.Pois(O.dataDT_rules, lambda = 3)
gstar1.N.Pois1.biyearly <- g.Pois(O.dataDT_rules, lambda = 1)
length(gstar1.N.Pois1.biyearly)
gstar2.N.p05 <- g.p(O.dataDT_rules, p = 0.5)
length(gstar2.N.p05)
O.dataDT_rules_Nstar <- cbind(O.dataDT_rules, gstar1.N.Pois3.yearly = gstar1.N.Pois3.yearly, gstar1.N.Pois1.biyearly = gstar1.N.Pois1.biyearly, gstar2.N.p05 = gstar2.N.p05)

# --------------------------------
# (IV) Estimate weights under observed (A,C,N)
# --------------------------------
# O.dataDT_rules[1:100,]
# O.dataDT_rules_Nstar[1:100, ]
gform.TRT <- "TI ~ CVD + highA1c + N.tminus1"
stratify.TRT <- list(
  TI=c("t == 0L",                                            # MODEL TI AT t=0
       "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
       "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
       "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
       # "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)"
      ))
gform.CENS <- c("C ~ highA1c")
stratify.CENS <- list(C=c("t < 16", "t == 16"))
gform.MONITOR <- "N ~ 1"
# stratify.TRT <- list(TI=c("t == 0L", "t > 0L & TI.tminus1 == 0L", "t > 0L & TI.tminus1 == 1L"))
# gform.CENS <- c("C + TI ~ highA1c + lastNat1", "N ~ highA1c + lastNat1 + C + TI")
# stratify.CENS <- list(C = NULL, TI = c("t == 0L", "t > 0"), N = c("t == 0L", "t > 0"))
# **** really want to define it like this ****
# gform.TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
#               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))
system.time(
  res <- stremr(data = O.dataDT_rules_Nstar, ID = "ID", t = "t",
                covars = c("highA1c", "lastNat1"),
                CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
                gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
                gform.TRT = gform.TRT, stratify.TRT = stratify.TRT,
                gform.MONITOR = gform.MONITOR,
                # gstar.TRT = "dlow",
                gstar.TRT = "dhigh",
                gstar.MONITOR = "gstar1.N.Pois3.yearly"
                # gstar.MONITOR = "gstar2.N.p05"
                )
  )
# Benchmark for N=50K:
#  user  system elapsed
# 9.269   2.121  11.345
names(res)

OData <- get_Odata(O.dataDT_rules_Nstar, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
modelfits.g0 <- get_fits(OData = OData, gform.CENS = gform.CENS, stratify.CENS = stratify.CENS, gform.TRT = gform.TRT, stratify.TRT = stratify.TRT, gform.MONITOR = gform.MONITOR)

wts.DT.dlow <- get_weights(modelfits.g0, OData, gstar.TRT = "dlow", gstar.MONITOR = "gstar1.N.Pois3.yearly")
survNP_ests.dlow <- get_survNP(wts.DT.dlow, OData)

wts.DT.dhigh <- get_weights(modelfits.g0, OData, gstar.TRT = "dhigh", gstar.MONITOR = "gstar1.N.Pois3.yearly")
survNP_ests.dhigh <- get_survNP(wts.DT.dhigh, OData)

OData$dat.sVar[]
# res$OData[1:100, ]
res[["IPW_estimates"]]
survS_dlow[13, "St"]-survS_dhigh[13, "St"]
# [1] 0.1655584 # CLOSE ENOUGH TO WHAT IS REPORTED on page 11, FIGURE 4: psi.N=0.173 (POISSON YEARLY)
survS_dhigh <- res[["IPW_estimates"]]
# new results dhigh:
#     t sum_Y_IPAW sum_all_IPAW          ht      m1ht        St
# 1   0   401.9918     49952.47 0.008047486 0.9919525 0.9919525
# 2   1   791.3740     50300.01 0.015733079 0.9842669 0.9763460
# 3   2  1564.1657     52158.19 0.029988877 0.9700111 0.9470665
# 4   3  2382.6140     56438.84 0.042215854 0.9577841 0.9070853
# 5   4  2530.2771     61485.51 0.041152411 0.9588476 0.8697566
# 6   5  2554.2661     66669.77 0.038312205 0.9616878 0.8364343
# 7   6  2614.1039     69971.42 0.037359593 0.9626404 0.8051854
# 8   7  2440.5948     72709.85 0.033566222 0.9664338 0.7781584
# 9   8  2050.3095     75970.86 0.026988105 0.9730119 0.7571574
# 10  9  1628.8897     74635.48 0.021824602 0.9781754 0.7406327
# 11 10  2109.0414     74228.24 0.028412926 0.9715871 0.7195892
# 12 11   564.9994     74996.10 0.007533718 0.9924663 0.7141680
# 13 12  1527.5176     75442.47 0.020247450 0.9797526 0.6997079
# 14 13   717.9634     77283.87 0.009289951 0.9907100 0.6932077
# 15 14   536.0347     80384.55 0.006668380 0.9933316 0.6885851
# 16 15  1929.0838     73191.46 0.026356679 0.9736433 0.6704363
# prev results:
# dhigh, gstar1.N.Pois3 (yearly):
#  1:  0   403.3415     50263.36 0.008024563 0.9919754 0.9919754
#  2:  1   846.9490     50265.85 0.016849393 0.9831506 0.9752613
#  3:  2  1581.6034     52216.13 0.030289557 0.9697104 0.9457210
#  4:  3  2387.5364     56128.93 0.042536644 0.9574634 0.9054932
#  5:  4  3087.5084     60770.58 0.050805971 0.9491940 0.8594888
#  6:  5  2890.4651     65220.91 0.044318071 0.9556819 0.8213979
#  7:  6  2872.1963     67695.78 0.042427997 0.9575720 0.7865476
#  8:  7  1660.7341     67806.23 0.024492353 0.9755076 0.7672832
#  9:  8  2433.1613     69533.91 0.034992443 0.9650076 0.7404341
# 10:  9  1947.7011     69949.34 0.027844455 0.9721555 0.7198171
# 11: 10   950.4954     69412.38 0.013693456 0.9863065 0.7099603
# 12: 11  1074.3250     69239.59 0.015516052 0.9844839 0.6989445
# 13: 12  1562.4440     68508.36 0.022806618 0.9771934 0.6830040
# 14: 13   764.1360     64659.01 0.011817935 0.9881821 0.6749323
# 15: 14   634.0271     64003.22 0.009906176 0.9900938 0.6682463
# 16: 15   210.8290     63797.74 0.003304647 0.9966954 0.6660380
survS_dlow <- res[["IPW_estimates"]]
# new results:
#     t sum_Y_IPAW sum_all_IPAW          ht      m1ht        St
# 1   0   410.4467     48624.36 0.008441174 0.9915588 0.9915588
# 2   1   496.1539     47567.03 0.010430625 0.9895694 0.9812162
# 3   2   487.9412     47098.47 0.010360022 0.9896400 0.9710508
# 4   3   357.0419     46909.43 0.007611304 0.9923887 0.9636599
# 5   4   397.0764     45988.38 0.008634276 0.9913657 0.9553394
# 6   5   640.4270     43429.09 0.014746500 0.9852535 0.9412514
# 7   6   718.1066     42472.89 0.016907411 0.9830926 0.9253373
# 8   7   419.0272     42755.80 0.009800477 0.9901995 0.9162686
# 9   8   595.5196     41638.50 0.014302136 0.9856979 0.9031640
# 10  9   414.2803     42053.07 0.009851368 0.9901486 0.8942666
# 11 10   842.6631     41863.64 0.020128761 0.9798712 0.8762661
# 12 11   228.3822     41177.78 0.005546248 0.9944538 0.8714061
# 13 12   287.3439     40782.04 0.007045845 0.9929542 0.8652663
# 14 13  1407.9491     40637.70 0.034646376 0.9653536 0.8352880
# 15 14   499.0048     40388.69 0.012355063 0.9876449 0.8249679
# 16 15  1396.3100     41049.10 0.034015607 0.9659844 0.7969062
# prev results:
#      t sum_Y_IPAW sum_all_IPAW          ht      m1ht        St
#  1:  0   427.1815     49609.03 0.008610962 0.9913890 0.9913890
#  2:  1   469.7380     49211.71 0.009545247 0.9904548 0.9819260
#  3:  2   602.5335     48892.29 0.012323692 0.9876763 0.9698250
#  4:  3   753.0877     48022.18 0.015682081 0.9843179 0.9546162
#  5:  4   972.1399     46747.80 0.020795415 0.9792046 0.9347645
#  6:  5   946.4126     44952.46 0.021053633 0.9789464 0.9150843
#  7:  6   918.9344     42631.18 0.021555454 0.9784445 0.8953593
#  8:  7   615.0077     42857.18 0.014350166 0.9856498 0.8825107
#  9:  8   875.2944     43870.59 0.019951738 0.9800483 0.8649031
# 10:  9   306.7152     42820.40 0.007162829 0.9928372 0.8587079
# 11: 10   996.5388     42677.21 0.023350607 0.9766494 0.8386566
# 12: 11   487.3777     41466.66 0.011753484 0.9882465 0.8287995
# 13: 12   892.2853     42468.93 0.021010308 0.9789897 0.8113861
# 14: 13   246.0274     41760.75 0.005891356 0.9941086 0.8066060
# 15: 14  1092.8948     40067.45 0.027276372 0.9727236 0.7846047
# 16: 15   392.7313     38375.61 0.010233879 0.9897661 0.7765751

# dlow, gstar2.N.p05 (no intervention on N.t)
# $St_ht_IPAW
#      t sum_Y_IPAW sum_all_IPAW          ht      m1ht        St
#  1:  0   428.6143     50040.85 0.008565289 0.9914347 0.9914347
#  2:  1   511.5621     49612.16 0.010311224 0.9896888 0.9812118
#  3:  2   592.3906     49100.69 0.012064813 0.9879352 0.9693737
#  4:  3   645.1490     48508.29 0.013299770 0.9867002 0.9564812
#  5:  4   773.8368     47862.66 0.016167861 0.9838321 0.9410170
#  6:  5   609.2330     47089.25 0.012937835 0.9870622 0.9288422
#  7:  6   605.9998     46479.86 0.013037900 0.9869621 0.9167321
#  8:  7   639.4173     45873.86 0.013938597 0.9860614 0.9039541
#  9:  8   602.7407     45233.98 0.013324953 0.9866750 0.8919090
# 10:  9   785.1512     44631.36 0.017591917 0.9824081 0.8762186
# 11: 10   507.8341     43846.12 0.011582189 0.9884178 0.8660701
# 12: 11   718.0278     43338.03 0.016568078 0.9834319 0.8517210
# 13: 12   769.8508     42619.86 0.018063194 0.9819368 0.8363362
# 14: 13   772.1379     41850.06 0.018450104 0.9815499 0.8209057
# 15: 14   693.4809     41077.71 0.016882168 0.9831178 0.8070470
# 16: 15   856.2705     40384.18 0.021203117 0.9787969 0.7899351

# dhigh, gstar2.N.p05 (no intervention on N.t)
#      t sum_Y_IPAW sum_all_IPAW          ht      m1ht        St
#  1:  0   408.1541     49999.60 0.008163148 0.9918369 0.9918369
#  2:  1  1050.3153     55609.10 0.018887472 0.9811125 0.9731036
#  3:  2  1682.9131     62829.03 0.026785597 0.9732144 0.9470384
#  4:  3  2456.9600     69437.27 0.035383878 0.9646161 0.9135285
#  5:  4  2797.7285     75717.63 0.036949500 0.9630505 0.8797741
#  6:  5  2817.0659     81031.54 0.034765053 0.9652349 0.8491887
#  7:  6  2400.5315     84203.29 0.028508762 0.9714912 0.8249794
#  8:  7  2273.0626     86276.06 0.026346388 0.9736536 0.8032442
#  9:  8  2305.1648     87122.22 0.026458977 0.9735410 0.7819911
# 10:  9  1656.7030     86812.55 0.019083680 0.9809163 0.7670679
# 11: 10  1543.2610     86149.69 0.017913714 0.9820863 0.7533268
# 12: 11  1425.5192     85064.46 0.016758105 0.9832419 0.7407025
# 13: 12  1259.8322     83913.38 0.015013485 0.9849865 0.7295820
# 14: 13  1160.6308     82824.52 0.014013131 0.9859869 0.7193582
# 15: 14  1278.4194     81767.37 0.015634836 0.9843652 0.7081112
# 16: 15  1177.7936     80517.99 0.014627708 0.9853723 0.6977532
res$OData
res$dataDT
summary(res$h_gN)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.006341 0.090590 0.322000 0.250300 0.340400 1.000000


# ------------------------------------------------------------------------------------------------------
# (IB) Simulate data with exacly the same variable ordering as in real data
# ------------------------------------------------------------------------------------------------------
DAGrm.realData <- set.DAG.realData()
DAGrm.realData <- set.DAG(DAGrm.realData)

# t1.wide <- system.time(O.data.w <- sim(DAG = DAGrm.realData, n = 10000, wide = TRUE))
# print(t1.wide)
#  user  system elapsed
# 3.339   0.102   3.445
# ncol(O.data.w)

t1.long <- system.time(O.data <- sim(DAG = DAGrm.realData, n = 10000, wide = FALSE))
print(t1.long)
# for 10K (with pre-alloc)
#  user  system elapsed
# 6.723   1.206   8.705

O.data[,"StudyID"] <- O.data[,"ID"]
O.data <- O.data[,!colnames(O.data)%in%c("ID","t")]
head(O.data[,1:5])
