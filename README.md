stremr
==========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/stremr)](http://cran.r-project.org/package=stremr)
[![](http://cranlogs.r-pkg.org/badges/stremr)](http://cran.rstudio.com/web/packages/stremr/index.html)
[![Travis-CI Build Status](https://travis-ci.org/osofr/stremr.svg?branch=master)](https://travis-ci.org/osofr/stremr)
[![Coverage Status](https://coveralls.io/repos/osofr/stremr/badge.svg?branch=master&service=github)](https://coveralls.io/github/osofr/stremr?branch=master)

The `stremr` R package implements the tools for streamlined analysis of  survival for treatment and monitoring events. In particular, it implements the Inverse Probability Weighted Estimator (IPW) of the hazard function in time-to-failure data. The user can specify static, dynamic or stochastic interventions on time-varying treatment, time-varying monitoring events. The input data needs to be in long format, with a specific **fixed** temporal ordering of the variables. The package supports modeling of right-censored data, where right-censoring variables can be coded as either binary or categorical (to represent all of the censoring events with one variable). Either representation allows for separate modeling of different right-censoring events. More than one variable can be used for coding each exposure and monitoring event and each of these can be binary, categorical or continuous. Separate models are fit for the observed censoring, exposure and monitoring mechanisms and each of these models can be fit by either pooling all subject-time data or by stratifying the model fits according to some arbitrary user-specified criteria. For example, the user may request that the exposure mechanism is modeled separately for each of the observed time-points (stratified by time). The same applies to modeling monitoring and censoring mechanisms. The output includes the estimated survival curve, which is obtained as a mapping from the estimated discrete hazard function. When several interventions for exposure/monitoring are specified, the package will produce one survival estimate for each intervention.

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

### Example with categorical censoring (3 levels)

Load the data:

```R
library("data.table")
library("magrittr")
data(OdataCatCENS)
OdataDT <- as.data.table(OdataCatCENS, key=c(ID, t))
# Indicator that the person has never been treated in the past:
OdataDT[, "barTIm1eq0" := as.integer(c(0, cumsum(TI)[-.N]) %in% 0), by = ID]
```


Regressions for modeling the exposure (TRT). Fit a separate model for TRT (stratify) for each of the following subsets:

```R
gform.TRT <- "TI ~ CVD + highA1c + N.tminus1"
stratify.TRT <- list(
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

Regressions for modeling the categorical censoring (CENS). Stratify the model fits by time-points (separate model for all t<16 and t=16):

```R
gform.CENS <- c("CatC ~ highA1c")
stratify.CENS <- list(CatC=c("t < 16", "t == 16"))
```

Regressions for modeling the monitoring regimen (MONITOR) and no stratification (pooling all observations over time):

```R
# Intercept only model, pooling across all time points t
gform.MONITOR <- "N ~ 1"
```

Define the probability of following a specific counterfactual monitoring regimen (the counterfactual probability of coming in for a visit is 0.1 at each time point):

```R
p <- 0.1 # probability of being monitored at each t is 0.1
OdataDT[, "gstar.N" := ifelse(N == 1L, eval(p), 1-eval(p))]

```

Define the indicator of following the counterfactual treatment rule (dynamic treatment rule) and perform estimation:

```R
# Define rule followers/non-followers for two rules: dlow & dhigh
res <- follow.rule.d.DT(OdataDT,
        theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
        CENS = "CatC", TRT = "TI", MONITOR = "N",
        rule.names = c("dlow", "dhigh")) %>%
# Merge rule definitions into main dataset:
  merge(OdataDT, ., by=c("ID", "t")) %>%
# Estimate hazard and survival for a rule "dhigh":
  stremr(gstar.TRT = "dhigh", gstar.MONITOR = "gstar.N",
        ID = "ID", t = "t", covars = c("highA1c", "lastNat1"),
        CENS = "CatC", gform.CENS = gform.CENS, stratify.CENS = stratify.CENS,
        TRT = "TI", gform.TRT = gform.TRT, stratify.TRT = stratify.TRT,
        MONITOR = "N", gform.MONITOR = gform.MONITOR, OUTCOME = "Y")

res$IPW_estimates
res$dataDT
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
