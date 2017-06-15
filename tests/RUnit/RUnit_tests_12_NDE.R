# # ---------------------------------------------------------------------------
# ## Test NDE argument "useonly_t_MONITOR":
# ## -- Will not intervene on g nodes that satisfy the logical expression provided to this arg
# # ---------------------------------------------------------------------------
# test.speedglm.stochastic.TMLE.NDE.1Kdata <- function() {
#   require("stremr")
#   options(stremr.verbose = FALSE)
#   options(width = 100)
#   `%+%` <- function(a, b) paste0(a, b)
#   require("data.table")
#   data(OdatDT_10K)
#   Odat_DT <- OdatDT_10K
#   Odat_DT <- Odat_DT[ID %in% (1:100), ]
#   # define intervention on N as 0101010101...
#   Odat_DT[, ("N.star.0101") := t%%2]
#   setkeyv(Odat_DT, cols = c("ID", "t"))

#   ## ---------------------------------------------------------------------------
#   ## Define some summaries (lags C[t-1], A[t-1], N[t-1])
#   ## ---------------------------------------------------------------------------
#   ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
#   lagnodes <- c("C", "TI", "N")
#   newVarnames <- lagnodes %+% ".tminus1"
#   Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
#   # indicator that the person has never been on treatment up to current t
#   Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
#   Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]
#   ## ----------------------------------------------------------------
#   ## IMPORT DATA
#   ## ----------------------------------------------------------------
#   set_all_stremr_options(estimator = "speedglm__glm")
#   OData <- importData(Odat_DT,
#                       ID = "ID", t = "t",
#                       covars = c("highA1c", "lastNat1", "lastNat1.factor"),
#                       CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome,
#                       remove_extra_rows = FALSE)
#   ## ------------------------------------------------------------------
#   ## Fit propensity scores for Treatment, Censoring & Monitoring
#   ## ------------------------------------------------------------------
#   gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
#   stratify_TRT <- list(
#     TI=c("t == 0L",                                            # MODEL TI AT t=0
#          "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
#          "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
#          "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
#         ))
#   gform_CENS <- c("C ~ highA1c + t")
#   gform_MONITOR <- "N ~ 1"

#   OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
#                           stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

#   ## ---------------------------------------------------------------------------------------------------------
#   ## IPW-KM with stochastic intervention on MONITOR
#   ## ---------------------------------------------------------------------------------------------------------
#   wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly")
#   surv1.stoch <- survNPMSM(wts.St.dlow, OData)
# #     time    sum_Y_IPW sum_all_IPAW     ht.NPMSM  St.NPMSM      ht.KM     St.KM           rule.name
# #  1:    0 3.136754e-01   22.1727277 1.414690e-02 0.9858531 0.06666667 0.9333333 gTI.dlowgPois3.yrly
# #  2:    1 0.000000e+00   12.5449959 0.000000e+00 0.9858531 0.00000000 0.9333333 gTI.dlowgPois3.yrly
# #  3:    2 0.000000e+00   18.0560124 0.000000e+00 0.9858531 0.00000000 0.9333333 gTI.dlowgPois3.yrly
# #  4:    3 0.000000e+00   11.3783485 0.000000e+00 0.9858531 0.00000000 0.9333333 gTI.dlowgPois3.yrly
# #  5:    4 0.000000e+00    1.8338345 0.000000e+00 0.9858531 0.00000000 0.9333333 gTI.dlowgPois3.yrly
# #  6:    5 7.486734e-03    0.9412119 7.954356e-03 0.9780113 0.07142857 0.8666667 gTI.dlowgPois3.yrly
# #  7:    6 0.000000e+00    0.9051833 0.000000e+00 0.9780113 0.00000000 0.8666667 gTI.dlowgPois3.yrly
# #  8:    7 5.949266e-01    1.0268491 5.793710e-01 0.4113799 0.07692308 0.8000000 gTI.dlowgPois3.yrly
# #  9:    8 0.000000e+00    0.5392287 0.000000e+00 0.4113799 0.00000000 0.8000000 gTI.dlowgPois3.yrly
# # 10:    9 0.000000e+00    0.5316701 0.000000e+00 0.4113799 0.00000000 0.8000000 gTI.dlowgPois3.yrly
# # 11:   10 0.000000e+00    0.3016605 0.000000e+00 0.4113799 0.00000000 0.8000000 gTI.dlowgPois3.yrly
# # 12:   11 0.000000e+00    0.4033243 0.000000e+00 0.4113799 0.00000000 0.8000000 gTI.dlowgPois3.yrly
# # 13:   12 2.084871e-03    0.6166494 3.380967e-03 0.4099890 0.08333333 0.7333333 gTI.dlowgPois3.yrly
# # 14:   13 0.000000e+00    0.5269700 0.000000e+00 0.4099890 0.00000000 0.7333333 gTI.dlowgPois3.yrly
# # 15:   14 4.731915e-08    0.5685754 8.322406e-08 0.4099890 0.09090909 0.6666667 gTI.dlowgPois3.yrly
# # 16:   15 0.000000e+00    0.7966520 0.000000e+00 0.4099890 0.00000000 0.6666667 gTI.dlowgPois3.yrly
# # 17:   16 0.000000e+00    0.0000000          NaN       NaN         NA        NA gTI.dlowgPois3.yrly

#   wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly")
#   surv2.stoch <- survNPMSM(wts.St.dhigh, OData)
# #     time  sum_Y_IPW sum_all_IPAW    ht.NPMSM  St.NPMSM      ht.KM     St.KM            rule.name
# #  1:    0 1.79840565    93.257054 0.019284393 0.9807156 0.01162791 0.9883721 gTI.dhighgPois3.yrly
# #  2:    1 0.52063384    82.780739 0.006289311 0.9745476 0.01250000 0.9760174 gTI.dhighgPois3.yrly
# #  3:    2 0.00000000    80.067556 0.000000000 0.9745476 0.00000000 0.9760174 gTI.dhighgPois3.yrly
# #  4:    3 3.04748751    67.052390 0.045449350 0.9302550 0.03389831 0.9429321 gTI.dhighgPois3.yrly
# #  5:    4 0.17038103    44.775942 0.003805192 0.9267152 0.06000000 0.8863562 gTI.dhighgPois3.yrly
# #  6:    5 3.47680389    41.222626 0.084342125 0.8485541 0.04878049 0.8431193 gTI.dhighgPois3.yrly
# #  7:    6 0.07590752    30.623653 0.002478722 0.8464508 0.05714286 0.7949410 gTI.dhighgPois3.yrly
# #  8:    7 4.93652926    37.463286 0.131769789 0.7349141 0.09375000 0.7204153 gTI.dhighgPois3.yrly
# #  9:    8 0.00000000    30.894441 0.000000000 0.7349141 0.00000000 0.7204153 gTI.dhighgPois3.yrly
# # 10:    9 4.85346504    25.464103 0.190600277 0.5948393 0.07142857 0.6689571 gTI.dhighgPois3.yrly
# # 11:   10 0.00000000    11.828032 0.000000000 0.5948393 0.00000000 0.6689571 gTI.dhighgPois3.yrly
# # 12:   11 0.18681931     8.504274 0.021967696 0.5817720 0.07692308 0.6174988 gTI.dhighgPois3.yrly
# # 13:   12 0.00000000     8.933052 0.000000000 0.5817720 0.00000000 0.6174988 gTI.dhighgPois3.yrly
# # 14:   13 0.00000000    10.928742 0.000000000 0.5817720 0.00000000 0.6174988 gTI.dhighgPois3.yrly
# # 15:   14 0.00000000     5.433945 0.000000000 0.5817720 0.00000000 0.6174988 gTI.dhighgPois3.yrly
# # 16:   15 0.00000000     8.716469 0.000000000 0.5817720 0.00000000 0.6174988 gTI.dhighgPois3.yrly
# # 17:   16 0.00000000     0.000000         NaN       NaN         NA        NA gTI.dhighgPois3.yrly

#   ## ---------------------------------------------------------------------------------------------------------
#   ## IPW-KM with static intervention on MONITOR
#   ## ---------------------------------------------------------------------------------------------------------
#   wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101")
#   surv1.stat <- survNPMSM(wts.St.dlow, OData)
# #     time sum_Y_IPW sum_all_IPAW   ht.NPMSM  St.NPMSM      ht.KM     St.KM           rule.name
# #  1:    0 0.2420811     16.80338 0.01440669 0.9855933 0.09090909 0.9090909 gTI.dlowN.star.0101
# #  2:    1 0.0000000     15.57885 0.00000000 0.9855933 0.00000000 0.9090909 gTI.dlowN.star.0101
# #  3:    2 0.0000000     13.39595 0.00000000 0.9855933 0.00000000 0.9090909 gTI.dlowN.star.0101
# #  4:    3 0.0000000     18.31408 0.00000000 0.9855933 0.00000000 0.9090909 gTI.dlowN.star.0101
# #  5:    4 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  6:    5 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  7:    6 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  8:    7 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  9:    8 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 10:    9 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 11:   10 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 12:   11 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 13:   12 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 14:   13 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 15:   14 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 16:   15 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 17:   16 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
#   wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101")
#   surv2.stat <- survNPMSM(wts.St.dhigh, OData)
# #     time sum_Y_IPW sum_all_IPAW   ht.NPMSM  St.NPMSM      ht.KM     St.KM            rule.name
# #  1:    0 0.9903319    49.190890 0.02013243 0.9798676 0.02222222 0.9777778 gTI.dhighN.star.0101
# #  2:    1 1.0892927    28.076497 0.03879732 0.9418513 0.04000000 0.9386667 gTI.dhighN.star.0101
# #  3:    2 0.0000000     8.452416 0.00000000 0.9418513 0.00000000 0.9386667 gTI.dhighN.star.0101
# #  4:    3 1.2218811     8.821175 0.13851682 0.8113891 0.16666667 0.7822222 gTI.dhighN.star.0101
# #  5:    4 0.0000000     6.935041 0.00000000 0.8113891 0.00000000 0.7822222 gTI.dhighN.star.0101
# #  6:    5 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# #  7:    6 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# #  8:    7 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# #  9:    8 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 10:    9 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 11:   10 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 12:   11 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 13:   12 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 14:   13 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 15:   14 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 16:   15 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 17:   16 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101

#   ## ---------------------------------------------------------------------------------------------------------
#   ## IPW-KM with static intervention on MONITOR under NDE assumption
#   ## -- Will not intervene on g nodes that satisfy the logical expression provided to this arg
#   ## ---------------------------------------------------------------------------------------------------------
#   wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
#   surv1.statNDE <- survNPMSM(wts.St.dlow, OData)
# #     time sum_Y_IPW sum_all_IPAW   ht.NPMSM  St.NPMSM      ht.KM     St.KM           rule.name
# #  1:    0 0.2420811     16.80338 0.01440669 0.9855933 0.09090909 0.9090909 gTI.dlowN.star.0101
# #  2:    1 0.0000000     15.57885 0.00000000 0.9855933 0.00000000 0.9090909 gTI.dlowN.star.0101
# #  3:    2 0.0000000     13.39595 0.00000000 0.9855933 0.00000000 0.9090909 gTI.dlowN.star.0101
# #  4:    3 0.0000000     18.31408 0.00000000 0.9855933 0.00000000 0.9090909 gTI.dlowN.star.0101
# #  5:    4 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  6:    5 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  7:    6 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  8:    7 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# #  9:    8 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 10:    9 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 11:   10 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 12:   11 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 13:   12 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 14:   13 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 15:   14 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 16:   15 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
# # 17:   16 0.0000000      0.00000        NaN       NaN         NA        NA gTI.dlowN.star.0101
#   wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
#   surv2.statNDE <- survNPMSM(wts.St.dhigh, OData)
# #     time sum_Y_IPW sum_all_IPAW   ht.NPMSM  St.NPMSM      ht.KM     St.KM            rule.name
# #  1:    0 0.9903319    49.190890 0.02013243 0.9798676 0.02222222 0.9777778 gTI.dhighN.star.0101
# #  2:    1 1.0892927    28.076497 0.03879732 0.9418513 0.04000000 0.9386667 gTI.dhighN.star.0101
# #  3:    2 0.0000000     8.452416 0.00000000 0.9418513 0.00000000 0.9386667 gTI.dhighN.star.0101
# #  4:    3 1.2218811     8.821175 0.13851682 0.8113891 0.16666667 0.7822222 gTI.dhighN.star.0101
# #  5:    4 0.0000000     6.935041 0.00000000 0.8113891 0.00000000 0.7822222 gTI.dhighN.star.0101
# #  6:    5 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# #  7:    6 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# #  8:    7 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# #  9:    8 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 10:    9 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 11:   10 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 12:   11 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 13:   12 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 14:   13 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 15:   14 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 16:   15 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101
# # 17:   16 0.0000000     0.000000        NaN       NaN         NA        NA gTI.dhighN.star.0101

#   ## ---------------------------------------------------------------------------------------------------------
#   ## TMLE / GCOMP with a stochastic intervention on MONITOR
#   ## ---------------------------------------------------------------------------------------------------------
#   t.surv <- c(0:2)
#   Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

#   gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
#   gcomp_est3$estimates[]
# #    est_name time  St.GCOMP St.TMLE ALLsuccessTMLE nFailedUpdates   type              IC.St fW_fit
# # 1:    GCOMP    0 0.9900000      NA          FALSE              1 pooled NA,NA,NA,NA,NA,NA,   NULL
# # 2:    GCOMP    1 0.9835187      NA          FALSE              2 pooled NA,NA,NA,NA,NA,NA,   NULL
# # 3:    GCOMP    2 0.9643009      NA          FALSE              3 pooled NA,NA,NA,NA,NA,NA,   NULL

#   # stratified modeling by rule followers only:
#   tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = TRUE)
#   tmle_est3$estimates[]
# #    est_name time St.GCOMP   St.TMLE ALLsuccessTMLE nFailedUpdates       type     SE.TMLE
# # 1:     TMLE    0       NA 0.9900000           TRUE              0 stratified 0.009949874
# # 2:     TMLE    1       NA 0.9805122           TRUE              0 stratified 0.015391437
# # 3:     TMLE    2       NA 0.9719458           TRUE              0 stratified 0.018072187

#   # pooling all observations (no stratification):
#   tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
#   tmle_est4$estimates[]
# #    est_name time St.GCOMP   St.TMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
# # 1:     TMLE    0       NA 0.9900000           TRUE              0 pooled 0.009949874
# # 2:     TMLE    1       NA 0.9805276           TRUE              0 pooled 0.015116604
# # 3:     TMLE    2       NA 0.9734881           TRUE              0 pooled 0.018452290

#   ## ---------------------------------------------------------------------------------------------------------
#   ## TMLE / GCOMP with a static intervention on MONITOR under NDE assumption
#   ## ---------------------------------------------------------------------------------------------------------
#   t.surv <- c(0:2)
#   Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

#   gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
#                             useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
#   gcomp_est3$estimates[]
# #    est_name time  St.GCOMP St.TMLE ALLsuccessTMLE nFailedUpdates   type              IC.St fW_fit
# # 1:    GCOMP    0 0.9900000      NA          FALSE              1 pooled NA,NA,NA,NA,NA,NA,   NULL
# # 2:    GCOMP    1 0.9679396      NA          FALSE              2 pooled NA,NA,NA,NA,NA,NA,   NULL
# # 3:    GCOMP    2 0.9616093      NA          FALSE              3 pooled NA,NA,NA,NA,NA,NA,   NULL

#   # stratified modeling by rule followers only:
#   tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
#                         useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
#   tmle_est3$estimates[]
# #    est_name time St.GCOMP   St.TMLE ALLsuccessTMLE nFailedUpdates       type     SE.TMLE
# # 1:     TMLE    0       NA 0.9900000           TRUE              0 stratified 0.009949874
# # 2:     TMLE    1       NA 0.9646167           TRUE              0 stratified 0.046448608
# # 3:     TMLE    2       NA 0.9434670           TRUE              0 stratified 0.062980775

#   # pooling all observations (no stratification):
#   tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
#                         useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
#   tmle_est4$estimates[]
# #    est_name time St.GCOMP   St.TMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
# # 1:     TMLE    0       NA 0.9900000           TRUE              0 pooled 0.009949874
# # 2:     TMLE    1       NA 0.9651137           TRUE              0 pooled 0.046793085
# # 3:     TMLE    2       NA 0.9423793           TRUE              0 pooled 0.064137543

#   ## ---------------------------------------------------------------------------------------------------------
#   ## TMLE / GCOMP with a stochastic intervention on MONITOR under NDE assumption
#   ## ---------------------------------------------------------------------------------------------------------
#   t.surv <- c(0:2)
#   Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
#   gcomp_est3 <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly",
#                             useonly_t_MONITOR = "gPois3.yrly == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
#   gcomp_est3$estimates[]
# #    est_name time  St.GCOMP St.TMLE ALLsuccessTMLE nFailedUpdates   type              IC.St fW_fit
# # 1:    GCOMP    0 0.9900000      NA          FALSE              1 pooled NA,NA,NA,NA,NA,NA,   NULL
# # 2:    GCOMP    1 0.9776050      NA          FALSE              2 pooled NA,NA,NA,NA,NA,NA,   NULL
# # 3:    GCOMP    2 0.9664007      NA          FALSE              3 pooled NA,NA,NA,NA,NA,NA,   NULL

#   # stratified modeling by rule followers only:
#   tmle_est3 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly",
#                         useonly_t_MONITOR = "gPois3.yrly == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
# #   tmle_est3$estimates[]
# # 1:     TMLE    0       NA 0.9900000           TRUE              0 stratified 0.009949874
# # 2:     TMLE    1       NA 0.9757197           TRUE              0 stratified 0.017747205
# # 3:     TMLE    2       NA 0.9900000           TRUE              0 stratified 0.013038014

#   # pooling all observations (no stratification):
#   tmle_est4 <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly",
#                         useonly_t_MONITOR = "gPois3.yrly == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
#   tmle_est4$estimates[]
# #    est_name time St.GCOMP   St.TMLE ALLsuccessTMLE nFailedUpdates   type     SE.TMLE
# # 1:     TMLE    0       NA 0.9900000           TRUE              0 pooled 0.009949874
# # 2:     TMLE    1       NA 0.9756022           TRUE              0 pooled 0.018267034
# # 3:     TMLE    2       NA 0.9787381           TRUE              0 pooled 0.018918870

# }
