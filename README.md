stremr
==========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/stremr)](http://cran.r-project.org/package=stremr)
[![](http://cranlogs.r-pkg.org/badges/stremr)](http://cran.rstudio.com/web/packages/stremr/index.html)
[![Travis-CI Build Status](https://travis-ci.org/osofr/stremr.svg?branch=master)](https://travis-ci.org/osofr/stremr)
[![Coverage Status](https://coveralls.io/repos/github/osofr/stremr/badge.svg?branch=master)](https://coveralls.io/github/osofr/stremr?branch=master)

The `stremr` R package implements the tools for streamlined analysis of  survival for treatment and monitoring events. In particular, it implements the Inverse Probability Weighted Estimator (IPW) of the hazard function in time-to-failure data. The user can specify static, dynamic or stochastic interventions on time-varying treatment, time-varying monitoring events. The input data needs to be in long format, with a specific **fixed** temporal ordering of the variables. The package supports modeling of right-censored data, where right-censoring variables can be coded as either binary or categorical (to represent all of the censoring events with one variable). Either representation allows for separate modeling of different right-censoring events. More than one variable can be used for coding each exposure and monitoring event and each of these can be binary, categorical or continuous. Separate models are fit for the observed censoring, exposure and monitoring mechanisms and each of these models can be fit by either pooling all subject-time data or by stratifying the model fitting according to some arbitrary user-specified criteria. For example, the user may request that the exposure mechanism is modeled separately for each of the observed time-points (stratified by time). The output includes the estimates of the discrete hazard function at each time-point along with the estimates of the discrete survival curve, obtained as a mapping from this hazard. When several interventions for exposure/monitoring are specified, the package will produce one survival estimate for each intervention.

* [Installing stremr](#Installation)
* [Documentation](#Documentation)
    * [Issue Tracking and Feature Requests](#IssueTracking) 
    * [List of Open Source Resources](#OpenSourceResources)
* [Example1](#Example1)
* [Creating Automatic Reports](#Reports)
* [Fitting Targeted Maximum Likelihood Estimation (TMLE) for longitudinal survival data](#TMLE)
* [Using H2O-3 Machine Learning Libraries](#H2OML)
* [Using Ensemble Learning (SuperLearner) based on H2O-3](#SuperLearner)

<a name="Installation"></a>
### Installation

<!-- To install the CRAN release version of `stremr`: 

```R
install.packages('stremr')
```
 -->

To install the development version (requires the `devtools` package):

```R
devtools::install_github('osofr/stremr', build_vignettes = FALSE)
```

<a name="Documentation"></a>
### Documentation

 For the general overview of the package:

```R
?stremr-package
```

For specific documentation on how to run `stremr()` function:
```R
?stremr
```

<!-- Once the package is installed, see the [vignette](http://cran.r-project.org/web/packages/stremr/vignettes/stremr_vignette.pdf), consult the internal package documentation and examples. 

* To see the vignette in R:

```R
vignette("stremr_vignette", package="stremr")
```

* To see all available package documentation:

```R
?stremr
help(package = 'stremr')
```

* To see the latest updates for the currently installed version of the package:

```r
news(package = "stremr")
```
 -->

<a name="Example1"></a>
### Example with categorical censoring (3 levels)

Load the data:

```R
library("data.table")
library("magrittr")
data(OdataCatCENS)
OdataDT <- as.data.table(OdataCatCENS, key=c(ID, t))
# Indicator that the person has never been treated in the past:
# ---------------------------------------------------------------------------
# Define some summaries (lags C[t-1], A[t-1], N[t-1])
# ---------------------------------------------------------------------------
ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
lagnodes <- c("C", "TI", "N")
newVarnames <- lagnodes %+% ".tminus1"
Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
# indicator that the person has never been on treatment up to current t
Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
```

Define two dynamic regimes (counterfactual treatment assignment under two rules (`dlow` & `dhigh`)):

```R
# Counterfactual TRT assignment for rule dlow (equivalent to always treated):
rule_name1 <- "dlow"
OdatDT[,"gTI." %+% rule_name1 := 1L]
# Counterfactual TRT assignment for dynamic rule dhigh -> start TRT only when I=1 (highA1c = 1)
rule_name2 <- "dhigh"
OdatDT_TIdhigh <- stremr::defineIntervedTRT(OdatDT, theta = 1, ID = ID,
                                          t = t, I = I, CENS = CENS, TRT = TRT, MONITOR = MONITOR,
                                          tsinceNis1 = "lastNat1",
                                          new.TRT.names = "gTI." %+% rule_name2)
OdatDT <- merge(OdatDT, OdatDT_TIdhigh, by=c(ID, t))
```


Define counterfactual monitoring probabilit(ies):

```R
# N^*(t) Bernoulli with P(N^*(t)=1)=p
g.p <- function(Odat, p) return(rep(p, nrow(Odat)))
# N^*(t) Poisson:
g.Pois <- function(Odat, lambda, lastNat1 = "lastNat1") {
  g.N <- 1 / (ppois(Odat[[lastNat1]], lambda, lower.tail = FALSE) / dpois(Odat[[lastNat1]], lambda) + 1)
  g.N[is.na(g.N)] <- 0
  return(g.N)
}
OdatDT <- OdatDT[, c("gPois3.yrly", "gPois3.biyrly", "gp05") := list(g.Pois(OdatDT, lambda = 3), g.Pois(OdatDT, lambda = 1), g.p(OdatDT, p = 0.5))][]
```

Regressions for modeling the exposure (TRT). Fit a separate model for TRT (stratify) for each of the following subsets:

```R
gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
stratify_TRT <- list(
  TI=c(
       # MODEL TI AT t=0
       "t == 0L",
       # MODEL TRT INITATION WHEN MONITORED
       "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",
       # MODEL TRT INITATION WHEN NOT MONITORED
       "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",
       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
       "(t > 0L) & (barTIm1eq0 == 1L)"
      ))
```

Regressions for modeling the categorical censoring (CENS). Stratify the model fits by time-points (separate model for all `t<16` and `t=16`):

```R
gform_CENS <- c("CatC ~ highA1c")
stratify_CENS <- list(CatC=c("t < 16", "t == 16"))
```

Regressions for modeling the monitoring regimen (MONITOR) and no stratification (pooling all observations over time):

```R
# Intercept only model, pooling across all time points t
gform_MONITOR <- "N ~ 1"
```

Define the probability of following a specific counterfactual monitoring regimen (the counterfactual probability of coming in for a visit is 0.1 at each time point):

```R
p <- 0.1 # probability of being monitored at each t is 0.1
OdataDT[, "gstar.N" := ifelse(N == 1L, eval(p), 1-eval(p))]

```

Define the indicator of counterfactual treatment (dynamic treatment rule): 

```R
# Define rule followers/non-followers for two rules: dlow & dhigh
res <- follow.rule.d.DT(OdataDT,
        theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
        CENS = "CatC", TRT = "TI", MONITOR = "N",
        rule.names = c("dlow", "dhigh")) %>%
# Merge rule definitions into main dataset:
  merge(OdataDT, ., by=c("ID", "t")) %>%
```

### Fit the propensity score models for censoring, exposure and monitoring:

```R
OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                        stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)
```



### Fit survival with non-parametric MSM  (IPTW-ADJUSTED KM):

```R
wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
survNPMSM(wts.St.dlow, OData)

wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
survNPMSM(wts.St.dhigh, OData)
```

### Turning the workflow into pipes:
```R
require("magrittr")
St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
            survNPMSM(OData) %$%
            IPW_estimates
St.dhigh
```

### Fitting IPW-adjusted MSM for hazard of the survival function

```R
MSM.IPAW <- survMSM(OData,
                    wts_data = list(dlow = wts.St.dlow, dhigh = wts.St.dhigh),
                    t_breaks = c(1:8,12,16)-1,
                    est_name = "IPAW", getSEs = FALSE)
MSM.IPAW
```


<a name="Reports"></a>
### Creating automatic reports:

html report:

```R
make_report_rmd(OData, MSM = MSM.IPAW, AddFUPtables = TRUE,
                RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = FALSE),
                WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                file.name = "sim.data.example")
```

pdf report:

```R
make_report_rmd(OData, MSM = MSM.IPAW, AddFUPtables = TRUE,
                RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = FALSE),
                WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                file.name = "sim.data.example", format = "pdf")

# omit extra modeling stuff (only coefficients):
make_report_rmd(OData, MSM = MSM.IPAW, RDtables = RDtables, file.path = report.path, only.coefs = TRUE, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
# skip modeling stuff alltogether:
make_report_rmd(OData, MSM = MSM.IPAW, RDtables = RDtables, file.path = report.path, skip.modelfits = TRUE, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
# skip RD tables by simply not including them:
make_report_rmd(OData, MSM = MSM.IPAW, file.path = report.path, skip.modelfits = TRUE, title = "Custom Report Title", author = "Oleg Sofrygin", y_legend = 0.95)
```

<a name="TMLE"></a>
### Fitting Targeted Maximum Likelihood Estimation (TMLE) for longitudinal survival data.

```R
t.surv <- c(10)
Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
params = list(fit.package = "speedglm", fit.algorithm = "GLM")
```

Stratified modeling by rule followers only:

```R
tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
tmle_est3
```

Pooling all observations (no stratification):

```R
tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
tmle_est4
```

RUN in PARALLEL seq-GCOMP & TMLE over t.surv (MUCH FASTER):

```R
require("doParallel")
registerDoParallel(cores = 40)
data.table::setthreads(1)
t.surv <- c(1,2,3,4,5,6,7,8,9,10)
Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

tmle_est_par1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE, parallel = TRUE)
tmle_est_par1
```

TMLE w/ h2o random forest:

```R
params = list(fit.package = "h2o", fit.algorithm = "RF", ntrees = 100,
              learn_rate = 0.05, sample_rate = 0.8,
              col_sample_rate = 0.8, balance_classes = TRUE)
t.surv <- c(10)
Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
tmle_est
```

### Citation

...
<!-- To cite `stremr` in publications, please use:
> Sofrygin O, van der Laan MJ, Neugebauer R (2015). *stremr: Simulating Longitudinal Data with Causal Inference Applications.* R package version 0.1.
 -->

### Funding

...
<!-- The development of this package was partially funded through internal operational funds provided by the Kaiser Permanente Center for Effectiveness & Safety Research (CESR). This work was also partially supported through a Patient-Centered Outcomes Research Institute (PCORI) Award (ME-1403-12506) and an NIH grant (R01 AI074345-07).
 -->

### Copyright
This software is distributed under the GPL-2 license.
