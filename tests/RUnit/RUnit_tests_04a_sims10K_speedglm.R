## --------------------------------------------------------------------------------------------------------
## Install data.table (most recent version)
# devtools::install_github('Rdatatable/data.table')
## --------------------------------------------------------------------------------------------------------
## Install stremr
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# ---------------------------------------------------------------------------

test.GCOMP.TMLE.10Kdata <- function() {
  `%+%` <- function(a, b) paste0(a, b)
  library("stremr")
  # options(stremr.verbose = TRUE)
  # options(gridisl.verbose = TRUE)
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)
  # set_all_stremr_options(estimator = "speedglm__glm")

  ## ---------------------------------------------------------------------------
  ## INSTALL CORRECT VERSIONS of data.table and stremr from github:
  ## ---------------------------------------------------------------------------
  # devtools::install_github('Rdatatable/data.table')
  # devtools::install_github('osofr/stremr', build_vignettes = FALSE)
  require("data.table")

  ## ---------------------------------------------------------------------------
  ## Test data set included in stremr:
  ## ---------------------------------------------------------------------------
  # head(O.data)
  data(OdatDT_10K)
  ID <- "ID"; t <- "t"; TRT <- "TI"; CENS <- "C"; MONITOR <- "N"; outcome <- "Y.tplus1"; I <- "highA1c";

  ## ---------------------------------------------------------------------------
  ## DEFINE SOME SUMMARIES (lags C[t-1], A[t-1], N[t-1])
  ## Might expand this in the future to allow defining arbitrary summaries
  ## ---------------------------------------------------------------------------
  # Odat_DT <- obsDTg05_1mil
  Odat_DT <- OdatDT_10K
  Odat_DT <- Odat_DT[ID %in% (1:100), ]
  lagnodes <- c("C", "TI", "N")
  newVarnames <- lagnodes %+% ".tminus1"
  Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
  # Indicator that the person has never been on treatment up to current t
  Odat_DT[, "barTIm1eq0" := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
  Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
  # Odat_DT[1:100, ]

  ## --------------------------------
  ## Define global options for stremr (which R packages to use for model fitting)
  ## --------------------------------
  # options(stremr.verbose = FALSE)
  # options(stremr.verbose = TRUE)

  # import data into stremr object:
  OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)

  ## --------------------------------
  ## Fitting the propensity scores for observed variables (A,C,N)
  ## --------------------------------
  # + N.tminus1
  gform_TRT <- "TI ~ CVD + highA1c"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))
  gform_CENS <- c("C ~ highA1c")
  stratify_CENS <- list(C=c("t < 16", "t == 16"))
  gform_MONITOR <- "N ~ 1"

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, stratify_CENS = stratify_CENS, gform_TRT = gform_TRT,
                                stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  # get IPW-adjusted and KM survival (with hazards over time)
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  St.dlow <- survNPMSM(wts.St.dlow, OData)
  St.dlow

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  St.dhigh <- survNPMSM(wts.St.dhigh, OData)
  St.dhigh

  ## ---------------------------------------------------------------------------------------------------------
  ## GCOMP AND TMLE w/ GLMs
  ## ---------------------------------------------------------------------------------------------------------
  # t.surv <- c(0,1,2,3,4,5,6,7,8,9,10)
  t.surv <- c(0:2)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  # stratified modeling by rule followers only:
  gcomp_est1 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE, stratify_by_last = FALSE)
  gcomp_est1$estimates[]
  # > gcomp_est1$estimates[]
  #    est_name  time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates       type
  #      <char> <num>     <num>  <lgcl>      <lgcl>         <lgcl>          <int>     <char>
  # 1:    GCOMP     1 0.9874747      NA          NA          FALSE              2 stratified
  # 2:    GCOMP     2 0.9710077      NA          NA          FALSE              3 stratified
  # 3:    GCOMP     3 0.9631855      NA          NA          FALSE              4 stratified
  # 4:    GCOMP    10 0.8660195      NA          NA          FALSE             11 stratified
  tmle_est1 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE, stratify_by_last = FALSE)
  tmle_est1$estimates[]
  #    est_name  time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates       type
  #      <char> <num>   <lgcl>     <num>      <lgcl>         <lgcl>          <int>     <char>
  # 1:     TMLE     1       NA 0.9884151          NA           TRUE              0 stratified
  # 2:     TMLE     2       NA 0.9662626          NA           TRUE              0 stratified
  # 3:     TMLE     3       NA 0.9578743          NA           TRUE              0 stratified
  # 4:     TMLE    10       NA 0.8613802          NA           TRUE              0 stratified

  gcomp_est1_last_t <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  gcomp_est1_last_t$estimates[]
  #    est_name  time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates       type
  #      <char> <num>     <num>  <lgcl>      <lgcl>         <lgcl>          <int>     <char>
  # 1:    GCOMP     1 0.9881351      NA          NA          FALSE              2 stratified
  # 2:    GCOMP     2 0.9775586      NA          NA          FALSE              3 stratified
  # 3:    GCOMP     3 0.9698591      NA          NA          FALSE              4 stratified
  # 4:    GCOMP    10 0.9028223      NA          NA          FALSE             11 stratified
  tmle_est1_last_t <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est1_last_t$estimates[]
  #    est_name  time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates       type
  #      <char> <num>   <lgcl>     <num>      <lgcl>         <lgcl>          <int>     <char>
  # 1:     TMLE     1       NA 0.9884149          NA           TRUE              0 stratified
  # 2:     TMLE     2       NA 0.9662519          NA           TRUE              0 stratified
  # 3:     TMLE     3       NA 0.9578642          NA           TRUE              0 stratified
  # 4:     TMLE    10       NA 0.8613870          NA           TRUE              0 stratified

  # pooling all observations (no stratification):
  gcomp_est2 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est2 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est2$estimates[]
  #    est_name  time  St.GCOMP St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type
  #      <char> <num>     <num>  <lgcl>      <lgcl>         <lgcl>          <int> <char>
  # 1:    GCOMP     1 0.9861232      NA          NA          FALSE              2 pooled
  # 2:    GCOMP     2 0.9751668      NA          NA          FALSE              3 pooled
  # 3:    GCOMP     3 0.9674501      NA          NA          FALSE              4 pooled
  # 4:    GCOMP    10 0.8964330      NA          NA          FALSE             11 pooled
  tmle_est2$estimates[]
  #    est_name  time St.GCOMP   St.TMLE St.iterTMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
  #      <char> <num>   <lgcl>     <num>      <lgcl>         <lgcl>          <int> <char>       <num>
  # 1:     TMLE     1       NA 0.9884138          NA           TRUE              0 pooled 0.001818848
  # 2:     TMLE     2       NA 0.9662455          NA           TRUE              0 pooled 0.004760497
  # 3:     TMLE     3       NA 0.9578574          NA           TRUE              0 pooled 0.005344252
  # 4:     TMLE    10       NA 0.8613797          NA           TRUE              0 pooled 0.009636681
  #                                                                             IC.St rule.name
  #                                                                            <list>    <char>
  # 1: -0.007387212,-0.006030586,-0.007387212, 0.003374826,-0.006030586,-0.006030586,  gTI.dlow
  # 2:       -0.01317465,-0.01027082,-0.01317465, 0.01702086,-0.01027082,-0.01027082,  gTI.dlow
  # 3:       -0.01573943,-0.01234550,-0.01573943, 0.03217811,-0.01234550,-0.01234550,  gTI.dlow
  # 4:       -0.03282305,-0.02735624,-0.03282305, 0.15442483,-0.02735624,-0.02735624,  gTI.dlow

  # stratified modeling by rule followers only:
  gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  gcomp_est3$estimates[]; tmle_est3$estimates[]

  # pooling all observations (no stratification):
  gcomp_est4 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est4$estimates[]; tmle_est4$estimates[]

  ## ------------------------------------------------------------------------
  ## RUN PARALLEL seq-GCOMP & TMLE over t.surv (MUCH FASTER)
  ## ------------------------------------------------------------------------
  # require("doParallel")
  # registerDoParallel(cores = 2)
  # data.table::setDTthreads(1)

  # gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # gcomp_est; tmle_est

  # gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  # gcomp_est; tmle_est
}

