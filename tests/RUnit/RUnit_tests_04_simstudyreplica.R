require("magrittr")
require("data.table")
require("stremr")

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

SimParams <- list(
  gcont       = list(name = "Continuous",         Nprob.t0 = 1,                                  Nprob.tplus = 1,                                     Nt.form = "N.t ~ 1",                  delta.shift = NULL),
  g05         = list(name = "Sporadic 0.5",       Nprob.t0 = 0.5,                                Nprob.tplus = 0.5,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  g035        = list(name = "Sporadic 0.35",      Nprob.t0 = 0.35,                               Nprob.tplus = 0.35,                                  Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  g02         = list(name = "Sporadic 0.2",       Nprob.t0 = 0.2,                                Nprob.tplus = 0.2,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  g01         = list(name = "Sporadic 0.1",       Nprob.t0 = 0.1,                                Nprob.tplus = 0.1,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
  # gSpoYear    = list(name = "Sporadic yearly",  Nprob.t0 = 0,                                  Nprob.tplus = substitute(N.sporadic.4(lastNat1[t])), Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = delta.shifts),
  # gSpoBiyear  = list(name = "Sporadic bi-yearly",Nprob.t0 = 0.25,                               Nprob.tplus = substitute(N.sporadic.2(lastNat1[t])), Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = 0.1),
  gPois3      = list(name = "Poisson yearly",     Nprob.t0 = substitute(N.Pois.3(lastNat1[t])),  Nprob.tplus = substitute(N.Pois.3(lastNat1[t])),     Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = delta.shifts, gradual.t.vec = 4),
  gPois1      = list(name = "Poisson bi-yearly",  Nprob.t0 = substitute(N.Pois.1(lastNat1[t])),  Nprob.tplus = substitute(N.Pois.1(lastNat1[t])),     Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = 0.1)
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
Nsize <- 10000
# set to TRUE to only run scenarios dealing with bias for g.N w and w/out past N(t), no intervention on N(t)
run_only_bias_noN <- FALSE
NoIntNt <- FALSE
# SimParams <- Sim.and.Est.params
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
O.dataDTrules <- add.gtheta.to.O(O.dataDT)
# --------------------------------
# (IV) Define monitoring regimen probabilit(ies):
# --------------------------------
gstar1.N.Pois3.yearly <- g.Pois(O.dataDTrules, lambda = 3)
gstar1.N.Pois1.biyearly <- g.Pois(O.dataDTrules, lambda = 1)
length(gstar1.N.Pois1.biyearly)
gstar2.N.p05 <- g.p(O.dataDTrules, p = 0.5)
length(gstar2.N.p05)

O.dataDTrules_Nstar <- O.dataDTrules[, gstar1.N.Pois3.yearly := gstar1.N.Pois3.yearly][, gstar1.N.Pois1.biyearly := gstar1.N.Pois1.biyearly][, gstar2.N.p05 := gstar2.N.p05]

# --------------------------------
# (IV) Estimate weights under observed (A,C,N)
# --------------------------------
# O.dataDTrules[1:100,]
# O.dataDTrules_Nstar[1:100, ]
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

OData <- get_Odata(O.dataDTrules_Nstar, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y")
OData <- get_fits(OData, gform.CENS = gform.CENS, stratify.CENS = stratify.CENS, gform.TRT = gform.TRT, stratify.TRT = stratify.TRT, gform.MONITOR = gform.MONITOR)

St.dlow <- get_weights(OData, gstar.TRT = "dlow", gstar.MONITOR = "gstar1.N.Pois3.yearly") %>%
           get_survNP(OData)  %$%
           IPW_estimates

St.dhigh <- get_weights(OData, gstar.TRT = "dhigh", gstar.MONITOR = "gstar1.N.Pois3.yearly") %>%
            get_survNP(OData) %$%
            IPW_estimates

OData$dat.sVar[]
# res$OData[1:100, ]
St.dlow[13, "St.IPTW"]-St.dhigh[13, "St.IPTW"]
St.list <- list(dlow = St.dlow[,"St.IPTW"], dhigh = St.dhigh[,"St.IPTW"])
# [1] 0.1779744 # CLOSE ENOUGH TO WHAT IS REPORTED on page 11, FIGURE 4: psi.N=0.173 (POISSON YEARLY)

wts.St.dlow <- get_weights(OData, gstar.TRT = "dlow", gstar.MONITOR = "gstar1.N.Pois3.yearly")
wts.St.dhigh <- get_weights(OData, gstar.TRT = "dhigh", gstar.MONITOR = "gstar1.N.Pois3.yearly")
wts.all.list <- list(dlow = wts.St.dlow, dhigh = wts.St.dhigh)

tjmin <- c(1:8,9,13)-1; tjmax <- c(1:8,12,16)-1

# MSM for hazard with regular weights:
MSM.IPAW <- get_survMSM(wts.all.list, OData, tjmin = tjmin, tjmax = tjmax, use.weights = TRUE, est.name = "IPAW")
# MSM for hazard with truncated weights:
MSM.trunc <- get_survMSM(wts.all.list, OData, tjmin = tjmin, tjmax = tjmax, use.weights = TRUE, trunc.weights = 20, est.name = "IPAWtrunc")
# crude MSM for hazard without any weights:
MSM.crude <- get_survMSM(wts.all.list, OData, tjmin = tjmin, tjmax = tjmax, use.weights = FALSE, est.name = "crude")


report.path <- "/Users/olegsofrygin/Dropbox/KP/monitoring_simstudy/stremr_test_report"
# report.path <- "/Users/olegsofrygin/Dropbox/KP/monitoring_simstudy/data_analysis/Reports"
# OData$emptydat.sVar
# save(list = c("OData", "MSM.IPAW", "MSM.trunc", "MSM.crude"), file = "MSM_results.RData")

# html doc:
make_report_rmd(OData, MSM = MSM.IPAW, file.path = report.path)
make_report_rmd(OData, MSM = MSM.trunc, file.path = report.path)
make_report_rmd(OData, MSM = MSM.crude, file.path = report.path)

# pdf doc:
make_report_rmd(OData, MSM = MSM.IPAW, format = "pdf", file.path = report.path)
make_report_rmd(OData, MSM = MSM.trunc, format = "pdf", file.path = report.path)
make_report_rmd(OData, MSM = MSM.crude, format = "pdf", file.path = report.path)
# word doc:
make_report_rmd(OData, MSM = MSM.IPAW, format = "word", file.path = report.path)

# make_report(OData, S.t = MSM.crude$St)

# system.time(
#   res <- stremr(data = O.dataDTrules_Nstar, ID = "ID", t = "t",
#                 covars = c("highA1c", "lastNat1"),
#                 CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y",
#                 gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
#                 gform.TRT = gform.TRT, stratify.TRT = stratify.TRT,
#                 gform.MONITOR = gform.MONITOR,
#                 # gstar.TRT = "dlow",
#                 gstar.TRT = "dhigh",
#                 gstar.MONITOR = "gstar1.N.Pois3.yearly"
#                 # gstar.MONITOR = "gstar2.N.p05"
#                 )
#   )
# # Benchmark for N=50K:
# #  user  system elapsed
# # 9.269   2.121  11.345
# names(res)

# # ------------------------------------------------------------------------------------------------------
# # (IB) Simulate data with exacly the same variable ordering as in real data
# # ------------------------------------------------------------------------------------------------------
# DAGrm.realData <- set.DAG.realData()
# DAGrm.realData <- set.DAG(DAGrm.realData)

# # t1.wide <- system.time(O.data.w <- sim(DAG = DAGrm.realData, n = 10000, wide = TRUE))
# # print(t1.wide)
# #  user  system elapsed
# # 3.339   0.102   3.445
# # ncol(O.data.w)

# t1.long <- system.time(O.data <- sim(DAG = DAGrm.realData, n = 10000, wide = FALSE))
# print(t1.long)
# # for 10K (with pre-alloc)
# #  user  system elapsed
# # 6.723   1.206   8.705

# O.data[,"StudyID"] <- O.data[,"ID"]
# O.data <- O.data[,!colnames(O.data)%in%c("ID","t")]
# head(O.data[,1:5])
