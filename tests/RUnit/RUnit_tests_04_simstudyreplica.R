# --------------------------------------------------------------------------------------------------------
# Install h2o (most recent version)
# --------------------------------------------------------------------------------------------------------
# if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
# if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
# # Next, we download packages that H2O depends on.
# if (! ("methods" %in% rownames(installed.packages()))) { install.packages("methods") }
# if (! ("statmod" %in% rownames(installed.packages()))) { install.packages("statmod") }
# if (! ("stats" %in% rownames(installed.packages()))) { install.packages("stats") }
# if (! ("graphics" %in% rownames(installed.packages()))) { install.packages("graphics") }
# if (! ("RCurl" %in% rownames(installed.packages()))) { install.packages("RCurl") }
# if (! ("jsonlite" %in% rownames(installed.packages()))) { install.packages("jsonlite") }
# if (! ("tools" %in% rownames(installed.packages()))) { install.packages("tools") }
# if (! ("utils" %in% rownames(installed.packages()))) { install.packages("utils") }
# # Now we download, install and initialize the H2O package for R.
# install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/rel-turchin/9/R")))

# ------------------------------------------------------------------------------------------------------
# SIMULATE, SAVE AND COMPRESS THE DATASET FROM THE EXAMPLE
# ------------------------------------------------------------------------------------------------------
notrun.save.example.data <- function() {
  require("data.table")
  require("simcausal")
  `%+%` <- function(a, b) paste0(a, b)
  # -----------------------------------------------------------
  # SIMULATION PARAMS:
  # -----------------------------------------------------------
  Nsize <- 1000000
  # Nsize <- 10000
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

  DAGrm <- set.DAG.rm.template()
  vecfun.add("N.Pois.1"); vecfun.add("N.Pois.3"); vecfun.add("N.sporadic.2"); vecfun.add("N.sporadic.4")
  #### Simulate with P(N(0)=1) = prob.t0 and P(N(t>0)=1) = prob.tplus
  if (!is.null(prob.t0)) {
    DAGrm <- DAGrm + node("N", t = 0, distr = "rbern", prob = .(prob.t0))
    DAGrm <- DAGrm + node("N", t = 1:16, distr = "rbern", prob = .(prob.tplus))
  }
  DAGrm <- set.DAG(DAGrm)
  Odat <- sim(DAG = DAGrm, n = Nsize, wide = FALSE, rndseed = 55466)
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
  OdatDT <- OdatDT[, c("gPois3.yrly", "gPois3.biyrly", "gp05") := list(g.Pois(OdatDT, lambda = 3), g.Pois(OdatDT, lambda = 1), g.p(OdatDT, p = 0.5))][]

  # --------------------------------
  # save data as csv
  # --------------------------------
  obsDTg05_1mil <- OdatDT
  data.table::fwrite(obsDTg05_1mil, "./obsDTg05_1mil.csv", turbo = TRUE, verbose = TRUE, na = "NA_h2o")

  # --------------------------------
  # save as compressed R file
  # --------------------------------
  # library("tools")
  # save(obsDTg05_1mil, compress = TRUE, file = "obsDTg05_1mil.rda", compression_level = 9)
  # resaveRdaFiles("./obsDTg05_1mil.rda", compress = "bzip2")
  # obsDTg05_10K <- OdatDT
  # save(obsDTg05_10K, compress = TRUE, file = "obsDTg05_10K.rda", compression_level = 9)
  # resaveRdaFiles("./obsDTg05_10K.rda", compress = "bzip2")
}

options(width = 100)
`%+%` <- function(a, b) paste0(a, b)
require("data.table")
obsDTg05_1mil <- data.table::fread(input = "./obsDTg05_1mil.csv", header = TRUE, na.strings = "NA_h2o")
# data(obsDTg05_1mil)
Odat_DT <- obsDTg05_1mil
setkeyv(Odat_DT, cols = c("ID", "t"))
head(Odat_DT)
nrow(Odat_DT)
# [1] 13945095

# ---------------------------------------------------------------------------
# DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
# ---------------------------------------------------------------------------
ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
lagnodes <- c("C", "TI", "N")
newVarnames <- lagnodes %+% ".tminus1"
Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
# indicator that the person has never been on treatment up to current t
Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

# ----------------------------------------------------------------
# IMPORT DATA
# ----------------------------------------------------------------
require("stremr")
options(stremr.verbose = TRUE)
# stremr_options(fit.package = "speedglm", fit.algorithm = "GLM")
stremr_options(fit.package = "h2o", fit.algorithm = "GLM")
# stremr_options(fit.package = "h2o", fit.algorithm = "RF")
# stremr_options(fit.package = "h2o", fit.algorithm = "GBM")

require("h2o")
h2o::h2o.init(nthreads = -1)
# h2o::h2o.shutdown(prompt = FALSE)

OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)

# --------------------------------
# Fit propensity scores for Treatment, Censoring & Monitoring
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
