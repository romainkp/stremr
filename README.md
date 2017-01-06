stremr
==========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/stremr)](https://CRAN.R-project.org/package=stremr)
[![](https://cranlogs.r-pkg.org/badges/stremr)](https://CRAN.R-project.org/package=stremr)
[![Travis-CI Build Status](https://travis-ci.org/osofr/stremr.svg?branch=master)](https://travis-ci.org/osofr/stremr)
[![Coverage Status](https://coveralls.io/repos/github/osofr/stremr/badge.svg)](https://coveralls.io/github/osofr/stremr)

Streamlined analysis of longitudinal time-to-event or time-to-failure data. Estimates the counterfactual discrete survival curve under static, dynamic and stochastic interventions on treatment (exposure) and monitoring events over time. Estimators (IPW, GCOMP, TMLE) adjust for *measured* time-varying confounding and informative right-censoring. Model fitting can be performed either with `glm` or [`H2O-3`](https://github.com/h2oai/h2o-3)machine learning libraries, including Ensemble Learning ([**SuperLearner**](https://github.com/h2oai/h2o-3/tree/master/h2o-r/ensemble)).

Currently implemented **estimators** include:
 - **Kaplan-Meier** Estimator. No adjustment for time-varying confounding or informative right-censoring.
 - [**Inverse Probability Weighted (IPW) Kaplan-Meier**](#survNPMSM). Also known as the Adjusted Kaplan Meier (AKME). Also known as the saturated (non-parametric) IPW-MSM estimator of the survival hazard. This estimator inverse weights each observation based on the exposure/censoring model fits (propensity scores).
 - [**Bounded Inverse Probability Weighted (B-IPW) Estimator of Survival**](#survDirectIPW). Estimates the survival directly (without hazard), also based on the exposure/censoring model fit (propensity scores).
 - [**Inverse Probability Weighted Marginal Structural Model (IPW-MSM)**](#survMSM) for the hazard function, mapped into survival. Currently only logistic regression is allowed where covariates are time-points and regime/rule indicators. This estimator is also based on the exposure/censoring model fit (propensity scores), but allows additional smoothing over multiple time-points and includes optional weight stabilization.
 - [**Sequential G-Computation (GCOMP)**](#GCOMPTMLE). Also known as the recursive G-Computation formula or Q-learning. Directly estimates the outcome model while adjusting for time-varying confounding. Estimation can be stratified by rule/regime followed or pooled across all rules/regimes.
 - [**Targeted Maximum Likelihood Estimator (TMLE)**](#GCOMPTMLE) for longitudinal data. Also known as the Targeted Minimum Loss-based Estimator. Doubly robust and semi-parametrically efficient estimator that targets the initial outcome model fits (GCOMP) with IPW.
 - **Iterative Targeted Maximum Likelihood Estimator (I-TMLE)** for longitudinal data. Fits sequential G-Computation and then iteratively performs targeting for all pooled Q's until convergence. 

**Input data**: 
 - Time-to-event (possibly) right-censored data has to be in long format.
 - Each row must contain a subject identifier (`ID`) and the integer indicator of the current time (`t`), e.g., day, week, month, year.
 - The package assumes that the temporal ordering of covariates in each row is **fixed** according to (`ID`, `t`, `L`,`C`,`A`,`N`,`Y`), where 
     * `L` -- Time-varying and baseline covariates.
     * `C` -- Indicators of right censoring events at time `t`; this can be either a single categorical or several binary columns.
     * `A` -- Exposure (treatment) at time `t`; this can be multivariate (more than one column) and each column can be binary, categorical or continuous.
     * `N` -- Indicator of being monitored at time point `t+1` (binary).
     * `Y` -- Time-to-event outcome (binary).
 - Note that the follow-up is assumed to end when either the outcome of interest (`Y[t]=1`) or right-censoring events are observed.
 - Categorical censoring can be useful for representing all of the censoring events with a single column (variable).

**Model fitting:**
 - Separate models are fit for the observed censoring, exposure and monitoring mechanisms.
 - Each model can be stratified (separate model is fit) by time or any other user-specified stratification criteria. Each strata is defined with by a single logical expression that selects specific observations/rows in the observed data (strata).
 -  By default, all models are fit using `GLM` with `binomial` family (logistic regression). 
 -  Alternatively, model fitting can be also performed with any machine learning algorithm implemented in `H2O-3` (faster distributed penalized `GLM`, `Random Forest`, `Gradient Boosting Machines` and `Deep Neural Network`).
 -  Finally, one can select the best model from an ensemble of H2O learners via cross-validation. Grid search (`h2o.grid`) allows for user-friendly model specification and fitting over multi-dimensional parameter space with various stopping criteria (random, discrete, max number of models, max time allocated, etc).
 -  The ensemble of many models can be combined into a single (more powerful) model with **SuperLearner** (`h2oEmsemble`). 

**Overview**:
* [Installing `stremr` and Documentation](#Installation)
* [Automated Reports](#Reports)
* [Example with Simulated Data](#Example1)
* [Sequential G-Computation (GCOMP) and Targeted Maximum Likelihood Estimation (TMLE) for longitudinal survival data](#GCOMPTMLE)
* [Machine Learning Algorithms](#H2OML)

<a name="Installation"></a>
### Installation and Documentation

<!-- To install the CRAN release version of `stremr`: 
```R
install.packages('stremr')
```
 -->

To install the development version (requires the `devtools` package):

```R
devtools::install_github('osofr/stremr', build_vignettes = FALSE)
```

For optimal performance, we also recommend installing the development version of `data.table`:
```R
devtools::install_github("Rdatatable/data.table")
```

For modeling with `H2O-3` machine learning libraries we recommend directly installing the latest version of the `h2o` R package ([can also see the instructions here](https://github.com/h2oai/h2o-3/tree/master/h2o-r#installation-from-within-r)):
```R
if ("package:h2o" %in% search()) detach("package:h2o", unload=TRUE)
if ("h2o" %in% rownames(installed.packages())) remove.packages("h2o")
# Next, download H2O package dependencies:
pkgs <- c("methods","statmod","stats","graphics","RCurl","jsonlite","tools","utils")
new.pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(new.pkgs)) install.packages(new.pkgs)
# Download and install the H2O package for R:
install.packages("h2o", type="source", repos=(c("https://s3.amazonaws.com/h2o-release/h2o/master/3636/R")))
```

Documentation with general overview of the package functions and datasets:
```R
?stremr-package
```

<!-- For specific documentation on how to run `stremr()` function:
```R
?stremr
```
 -->

To obtain documentation for specific relevant functions in `stremr` package:
```R
?importData
?fitPropensity
?getIPWeights
?survDirectIPW
?survNPMSM
?survMSM
?fitSeqGcomp
?fitTMLE
?fitIterTMLE
```

<a name="Reports"></a>
### Automated Reports:

The following is an example of a function call that produces an automated `html` report shown below. For a pdf report just set the argument `format = "pdf"`.
```R
  make_report_rmd(OData, NPMSM = list(surv1, surv2), 
                  MSM = MSM.IPAW, 
                  GCOMP = list(gcomp_est1, gcomp_est2), 
                  TMLE = list(tmle_est_par1, tmle_est_par2),
                  AddFUPtables = TRUE, RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Author Name", y_legend = 0.99, x_legend = 9.5)
```

![gif](https://cloud.githubusercontent.com/assets/6721358/18609476/d9b4db74-7cb7-11e6-9ca6-aacf0b70ca4c.gif)


<a name="Example1"></a>
### Example with Simulated Data

Load the data:

```R
require("stremr")
require("data.table")
data(OdataNoCENS)
OdataDT <- as.data.table(OdataNoCENS, key=c(ID, t))
```

Define some summaries (lags):
```R
OdataDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
OdataDT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
```

Import input data into `stremr` object `DataStorageClass` and define relevant covariates:
```R
OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "N.tminus1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
```

Define counterfactual exposures. In this example we define one intervention as always treated  and another as never treated. Such intervention can be defined conditionally on other variables (dynamic intervention). Similarly, one can define the intervention as a probability that the counterfactual exposure is 1 at each time-point `t` (for stochastic interventions).
```R
OdataDT[, ("TI.set1") := 1L]
OdataDT[, ("TI.set0") := 0L]
```

Regressions for modeling the propensity scores for censoring (`CENS`), exposure (`TRT`) and monitoring (`MONITOR`). By default, each of these propensity scores is fit with a common model that pools across all available time points (smoothing over time).
```R
gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
gform_MONITOR <- "N ~ 1"
```

Stratification, that is, fitting separate models for different time-points, is enabled with logical expressions in arguments `stratify_...` (see `?fitPropensity`). For example, the logical expression below states that we want to fit the censoring mechanism with a separate model for time point 16, while pooling with a common model fit over time-points 0 to 15. Any logical expression can be used to define such stratified modeling. This can be similarly applied to modeling the exposure mechanism (`stratify_TRT`) and the monitoring mechanism (`stratify_MONITOR`).
```R
stratify_CENS <- list(C=c("t < 16", "t == 16"))
```

Fit the propensity scores for censoring, exposure and monitoring:
```R
OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR, stratify_CENS = stratify_CENS)
```

<a name="survNPMSM"></a>Estimate survival based on non-parametric MSM (IPTW-ADJUSTED KM):
```R
require("magrittr")
AKME.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
             survNPMSM(OData) %$%
             IPW_estimates
AKME.St.1
```

<a name="survDirectIPW"></a>Estimate survival with bounded IPW:
```R
IPW.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
            survDirectIPW(OData)
IPW.St.1[]
```

<a name="survMSM"></a>Estimate hazard with IPW-MSM then map into survival estimate. Using two regimens and smoothing over two intervals of time-points:
```R
wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, t_breaks = c(1:8,12,16)-1,)
survMSM_res$St
```

<a name="GCOMPTMLE"></a>
### Sequential G-Computation (GCOMP) and Targeted Maximum Likelihood Estimation (TMLE) for longitudinal time-to-event data.

Define time-points of interest, regression formulas and software to be used for fitting the sequential outcome models:
```R
t.surv <- c(0:15)
Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
params = list(fit.package = "speedglm", fit.algorithm = "glm")
```

G-Computation (pooled):
```R
gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
```

Targeted Maximum Likelihood Estimation (TMLE) (stratified):
```R
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE)
tmle_est[]
```

To parallelize estimation over several time-points (`t.surv`) for either GCOMP or TMLE use argument `parallel = TRUE`:
```R
require("doParallel")
registerDoParallel(cores = 40)
data.table::setthreads(1)
tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = TRUE, parallel = TRUE)
```

<a name="H2OML"></a>
### Machine Learning Algorithms

To perform all modeling with `H2O-3` distributed Random Forest algorithm just set the global package options `fit.package = "h2o"` and `fit.algorithm = "randomForest"` prior to calling any fitting function:
```R
set_all_stremr_options(fit.package = "h2o", fit.algorithm = "randomForest")

require("h2o")
h2o::h2o.init(nthreads = -1)

OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR, stratify_CENS = stratify_CENS)
```

Other available algorithms are `H2O-3` Gradient Boosting Machines (`fit.algorithm = "gbm"`), distributed GLM (including LASSO and Ridge) (`fit.algorithm = "glm"`) and Deep Neural Nets (`fit.algorithm = "deeplearning"`).

Use arguments `params_...` in `fitPropensity()` and `params_Q` in `fitSeqGcomp()` and `fitTMLE()` to pass various tuning parameters and select different algorithms for different models:
```R
params_TRT = list(fit.package = "h2o", fit.algorithm = "gbm", ntrees = 50, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)
params_CENS = list(fit.package = "speedglm", fit.algorithm = "glm")
params_MONITOR = list(fit.package = "speedglm", fit.algorithm = "glm")
OData <- fitPropensity(OData,
          gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, params_CENS = params_CENS,
          gform_TRT = gform_TRT, params_TRT = params_TRT,
          gform_MONITOR = gform_MONITOR, params_MONITOR = params_MONITOR)
```

Running TMLE based on the previous fit of the propensity scores. Also applying Random Forest to estimate the sequential outcome model:
```R
params_Q = list(fit.package = "h2o", fit.algorithm = "randomForest", ntrees = 100, learn_rate = 0.05, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)

tmle_est <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, params_Q = params_Q, stratifyQ_by_rule = TRUE)
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
