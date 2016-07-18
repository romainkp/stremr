require("data.table")
`%+%` <- function(a, b) paste0(a, b)

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
    # node("TI", t = 0, distr = "rbern", prob = 1) +
    node("TI", t = 0, distr = "rbern", prob = ifelse(highA1c[0] == 0, ifelse(CVD[0] == 1, 0.5, 0.1), ifelse(CVD[0] == 1, 0.9, 0.5))) +
    node("C", t = 0, distr = "rbern", prob = 0, EFU = TRUE) +  # no censoring in first bin
    node("N", t = 0, distr = "rbern", prob = 1) +
    node("Y", t = 1:16, distr = "rbern", prob = plogis(-6.5 + 1*CVD[0] + 4*highA1c.UN[t-1] + 0.05*timelowA1c.UN[t-1]),  EFU = TRUE) +  # Z(t)
    node("lastNat1", t = 1:16, distr = "rconst", const = ifelse(N[t-1] == 0, lastNat1[t-1] + 1, 0)) +  # Z(1)  just a function of past \bar{N}(t-1) - 0 probs current N at 1,  1 probs previous N a 1,  2 probs  the one before the previous was at 1,  etc.
    node("highA1c.UN", t = 1:16, distr = "rbern", prob = ifelse(TI[t-1] == 1, 0.1, ifelse(highA1c.UN[t-1] == 1, 0.9, min(1, 0.1 + t/16)))) +  # I(t)
    node("highA1c", t = 1:16, distr = "rbern", prob = ifelse(N[t-1] == 1, highA1c.UN[t], highA1c[t-1])) + # I*(m)=I(m)*T(m-1)  (I actually replace I*(m)=0 with when T(m-1)=0 with last value carried forward,  i.e. I*(m-1)
    node("timelowA1c.UN", t=1:16, distr="rnorm", mean=sum(1-highA1c.UN[0:t]),  sd=0) +  # Z(m)
    # node("TI", t = 1:16, distr = "rbern", prob = 1) +
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

# --------------------------------
# Define event indicator and right-censoring indicator, define total follow-up length
# --------------------------------
define_indicators <- function(O.data, ID = "ID", t = "t", TRT = "TI", CENS = "C", MONITOR = "N", I = "highA1c") {
  # Add indicator Delta=I(T<=C) (Y is 1 at some point for the subject):
  O.dataDT <- data.table(O.data, key = c(ID, t))
  EVENT_IND <- "Delta";
  O.dataDT[,(EVENT_IND) := as.integer(any(get(outcome) %in% 1)), by=eval(ID)]
  # Add indicator Delta=I(T>C) (subject was right-censored at some point):
  CENS_IND <- "AnyCensored"; noCENScat <- 0L; CENS <- c("C")
  O.dataDT[, (CENS_IND) := FALSE, by=eval(ID)]
  for (Cvar in CENS) {
    O.dataDT[, (CENS_IND) := get(CENS_IND) | any(!get(Cvar) %in% c(eval(noCENScat),NA)), by = eval(ID)]
  }
  O.dataDT[, (CENS_IND) := as.integer(get(CENS_IND))]
  ## Add variable that indicates if TI initiation occured previously, relevant for model of TI continuation
  TIcovarname <- "barTIm1eq0"
  O.dataDT[, (TIcovarname) := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  # O.dataDT[1:100,]
  return(O.dataDT)
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

  # counterfactual TRT for rule dlow (equivalent to always treated):
  O.dataDT[,"TI.gstar." %+% rule_name1 := 1L]

  # Add indicator for rule (dhigh) -> start TRT only when I=1 (highA1c)
  rule_name2 <- "dhigh"
  O.dataDT_dhigh <- stremr::defineTRTrules(O.dataDT, theta = 1, ID = ID,
                                              t = t, I = I, CENS = CENS, TRT = TRT,
                                              MONITOR = MONITOR,
                                              tsinceNis1 = "lastNat1",
                                              rule.names = rule_name2)
  O.dataDT_dhigh[, (rule_name2) := as.integer(get(rule_name2))]
  O.dataDT <- merge(O.dataDT, O.dataDT_dhigh, by=c(ID, t))

  O.dataDT_TIdhigh <- stremr::defineIntervedTRT(O.dataDT, theta = 1, ID = ID,
                                            t = t, I = I, CENS = CENS, TRT = TRT,
                                            MONITOR = MONITOR,
                                            tsinceNis1 = "lastNat1",
                                            new.TRT.names = "TI.gstar." %+% rule_name2)
  O.dataDT <- merge(O.dataDT, O.dataDT_TIdhigh, by=c(ID, t))
  O.dataDT[1:100, ]

  return(O.dataDT)
}


# ------------------------------------------------------------------------------------------------------
# SIMULATE, SAVE AND COMPRESS THE DATASET FROM THE EXAMPLE
# ------------------------------------------------------------------------------------------------------
notrun.save.example.data <- function() {
  require("simcausal")
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
  # -----------------------------------------------------------
  # SIMULATION PARAMS:
  # -----------------------------------------------------------
  # Nsize <- 200000
  Nsize <- 100000
  # Nsize <- 10000
  # set to TRUE to only run scenarios dealing with bias for g.N w and w/out past N(t), no intervention on N(t)
  # SimParams <- Sim.and.Est.params
  # override the scenario name to run for just one scenario:
  scen.names <- "g05"
  delta.shifts <- c(seq(-0.9, -0.1, by=0.1), seq(0.1, 0.9, by=0.1))
  # scen.names <- c("gPois3", "gPois1")
  # scen.names <- c("g05", "gPois3")
  # scen.names <- c("gSpoYear", "gPois3")
  SimParams <- list(
    gcont       = list(name = "Continuous",         Nprob.t0 = 1,                                  Nprob.tplus = 1,                                     Nt.form = "N.t ~ 1",                  delta.shift = NULL),
    g05         = list(name = "Sporadic 0.5",       Nprob.t0 = 0.5,                                Nprob.tplus = 0.5,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
    g035        = list(name = "Sporadic 0.35",      Nprob.t0 = 0.35,                               Nprob.tplus = 0.35,                                  Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
    g02         = list(name = "Sporadic 0.2",       Nprob.t0 = 0.2,                                Nprob.tplus = 0.2,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
    g01         = list(name = "Sporadic 0.1",       Nprob.t0 = 0.1,                                Nprob.tplus = 0.1,                                   Nt.form = "N.t ~ 1",                  delta.shift = 0.1),
    gPois3      = list(name = "Poisson yearly",     Nprob.t0 = substitute(N.Pois.3(lastNat1[t])),  Nprob.tplus = substitute(N.Pois.3(lastNat1[t])),     Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = delta.shifts, gradual.t.vec = 4),
    gPois1      = list(name = "Poisson bi-yearly",  Nprob.t0 = substitute(N.Pois.1(lastNat1[t])),  Nprob.tplus = substitute(N.Pois.1(lastNat1[t])),     Nt.form = "N.t ~ as.factor(lastN.t)", delta.shift = 0.1)
  )
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
  # O.data.simstudy <- sim(DAG = DAGrm, n = Nsize, wide = TRUE, LTCF = "Y", rndseed = 55466)
  # head(O.data.simstudy)
  # colMeans(O.data.simstudy[,"Y_"%+%c(1:16)])

  O.data.simstudy[O.data.simstudy[,"t"]%in%16,"lastNat1"] <- NA
  O.data.simstudy.g05 <- O.data.simstudy
  # O.data.simstudy.g05 <- O.data.simstudy[,!names(O.data.simstudy)%in%c("highA1c.UN", "timelowA1c.UN")]

  save(O.data.simstudy.g05, compress = TRUE, file = "O.data.simstudy.g05.rda", compression_level = 9)
  library("tools")
  resaveRdaFiles("./O.data.simstudy.g05.rda", compress = "bzip2")
}

data(O.data.simstudy.g05)
O.data <- O.data.simstudy.g05
head(O.data)
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
# get the likelihood under N(t) Bernoulli with P(N(t)=1)=p:
g.p <- function(O.dataDT, p, ID = "ID", OUTCOME = "Y", lastNat1 = "lastNat1", MONITOR = "N"){
  g.N <- rep(p, nrow(O.dataDT))
  # g.N <- ifelse(O.dataDT[[MONITOR]] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
  return(g.N)
}
# get the likelihood under N(t) Poisson:
g.Pois <- function(O.dataDT, lambda, ID = "ID", OUTCOME = "Y", lastNat1 = "lastNat1", MONITOR = "N"){
  g.N <- 1 / (ppois(O.dataDT[[lastNat1]], lambda, lower.tail = FALSE) / dpois(O.dataDT[[lastNat1]], lambda) + 1)
  g.N[is.na(g.N)] <- 0
  # g.N <- ifelse(O.dataDT[[MONITOR]] == 1, g.N, 1 - g.N) # will make it NA for rows when Y.t=1
  return(g.N)
}

gstar1.N.Pois3.yearly <- g.Pois(O.dataDTrules, lambda = 3)
gstar1.N.Pois1.biyearly <- g.Pois(O.dataDTrules, lambda = 1)
gstar2.N.p05 <- g.p(O.dataDTrules, p = 0.5)
O.dataDTrules_Nstar <- O.dataDTrules[, gstar1.N.Pois3.yearly := NULL][, gstar1.N.Pois1.biyearly := NULL][, gstar2.N.p05 := NULL]
O.dataDTrules_Nstar <- O.dataDTrules[, gstar1.N.Pois3.yearly := gstar1.N.Pois3.yearly][, gstar1.N.Pois1.biyearly := gstar1.N.Pois1.biyearly][, gstar2.N.p05 := gstar2.N.p05]
# O.dataDTrules_Nstar <- O.dataDTrules
# ---------------------------------------------------------------------------
# DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
# Might expand this in the future to allow defining arbitrary summaries
# ---------------------------------------------------------------------------
lagnodes <- c("C", "TI", "N")
newVarnames <- lagnodes %+% ".tminus1"
O.dataDTrules_Nstar[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
# -------------------------------------------------------------------------------------------
# Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
# -------------------------------------------------------------------------------------------
OUTCOME <- "Y"
shifted.OUTCOME <- OUTCOME%+%".tplus1"
O.dataDTrules_Nstar[, (shifted.OUTCOME) := shift(get(OUTCOME), n = 1L, type = "lead"), by = ID]
O.dataDTrules_Nstar <- O.dataDTrules_Nstar[!get(OUTCOME)%in%1,]
# O.dataDTrules_Nstar[1:50]

# --------------------------------
# (IV) Estimate weights under observed (A,C,N)
# --------------------------------
gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
stratify_TRT <- list(
  TI=c("t == 0L",                                            # MODEL TI AT t=0
       "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
       "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
       "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
      ))
gform_CENS <- c("C ~ highA1c")
stratify_CENS <- list(C=c("t < 16", "t == 16"))
gform_MONITOR <- "N ~ 1"
# **** really want to define it like this ****
# gform_TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
#               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))

# ----------------------------------------------------------------
# TESTING IPW and IPW-MSM. Fitting with speedglm, h2oglm, RF and GBM
# ----------------------------------------------------------------
options(stremr.verbose = TRUE)
require("h2o")
# h2o::h2o.init(nthreads = 2)
h2o::h2o.init()
stremr_options(fit.package = "speedglm", fit.algorithm = "GLM")
# stremr_options(fit.package = "h2o", fit.algorithm = "GLM"); model <- "h2o.GLM"
# stremr_options(fit.package = "h2o", fit.algorithm = "RF"); model <- "h2o.RF"
# stremr_options(fit.package = "h2o", fit.algorithm = "GBM"); model <- "h2o.GBM"
# h2o::h2o.shutdown(prompt = FALSE)

OData <- importData(O.dataDTrules_Nstar, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = shifted.OUTCOME)
# OData$fast.load.to.H2O()
# OData$H2O.dat.sVar
OData$dat.sVar

params_CENS = list()
params_TRT = list()
params_MONITOR = list()

OData <- fitPropensity(OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT,
                              stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
                              params_CENS = params_CENS, params_TRT = params_TRT, params_MONITOR = params_MONITOR)

require("magrittr")
# OData$dat.sVar
St.dlow <- getIPWeights(OData, gstar_TRT = "TI.gstar.dlow", gstar_MONITOR = "gstar1.N.Pois3.yearly") %>%
           survNPMSM(OData)  %$%
           IPW_estimates
St.dlow

St.dhigh <- getIPWeights(OData, gstar_TRT = "TI.gstar.dhigh", gstar_MONITOR = "gstar1.N.Pois3.yearly") %>%
            survNPMSM(OData) %$%
            IPW_estimates
St.dhigh

St.dlow[13, "St.IPTW"]-St.dhigh[13, "St.IPTW"] # [1] 0.1260063
St.list <- list(dlow = St.dlow[,"St.IPTW"], dhigh = St.dhigh[,"St.IPTW"])

wts.St.dlow <- getIPWeights(OData, gstar_TRT = "TI.gstar.dlow")
# wts.St.dlow <- getIPWeights(OData, gstar_TRT = "dlow")
wts.St.dhigh <- getIPWeights(OData, gstar_TRT = "TI.gstar.dhigh")
# wts.St.dhigh <- getIPWeights(OData, gstar_TRT = "dhigh")
# wts.St.dlow <- getIPWeights(OData, gstar_TRT = "dlow", gstar_MONITOR = "gstar1.N.Pois3.yearly")
# wts.St.dhigh <- getIPWeights(OData, gstar_TRT = "dhigh", gstar_MONITOR = "gstar1.N.Pois3.yearly")
wts.all.list <- list(dlow = wts.St.dlow, dhigh = wts.St.dhigh)
wts.all <- rbindlist(wts.all.list)

tjmin <- c(1:8,9,13)-1; tjmax <- c(1:8,12,16)-1
# MSM for hazard with regular weights:
MSM.IPAW <- survMSM(OData, wts_data = wts.all, t_breaks = tjmax,
                    # tjmin = tjmin, tjmax = tjmax,
                    use_weights = TRUE, est_name = "IPAW", getSEs = TRUE)
report.path <- "/Users/olegsofrygin/Dropbox/KP/monitoring_simstudy/stremr_test_report"

# print everything:
# make_report_rmd(OData, MSM = MSM.IPAW, AddFUPtables = TRUE,
#                 RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
#                 WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
#                 file.name = model%+%"_sim.data", file.path = report.path, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)

# make_report_rmd(OData, MSM = MSM.IPAW,
#                 RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
#                 WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
#                 file.name = model%+%"_sim.data", file.path = report.path, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)

# # omit extra h2o model output stuff (only coefficients):
# make_report_rmd(OData, MSM = MSM.IPAW, RDtables = RDtables, file.path = report.path, only.coefs = TRUE, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
# # skip modeling stuff alltogether:
# make_report_rmd(OData, MSM = MSM.IPAW, RDtables = RDtables, file.path = report.path, skip.modelfits = TRUE, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
# # skip RD tables by simply not including them:
# make_report_rmd(OData, MSM = MSM.IPAW, file.path = report.path, skip.modelfits = TRUE, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)

# make_report_rmd(OData, MSM = MSM.IPAW, RDtables = RDtables, format = "pdf", file.path = report.path, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
# make_report_rmd(OData, MSM = MSM.IPAW, RDtables = RDtables, format = "word",file.path = report.path, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
