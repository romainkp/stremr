## --------------------------------------------------------------------------------------------------------
## Install data.table (most recent version)
# devtools::install_github('Rdatatable/data.table')
## --------------------------------------------------------------------------------------------------------
## Install stremr
# devtools::install_github('osofr/stremr', build_vignettes = FALSE)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test speedglm
# ---------------------------------------------------------------------------
test.speedglm.allestimators10Kdata <- function() {
  # options(stremr.verbose = FALSE)
  options(stremr.verbose = TRUE)
  options(width = 100)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  # obsDTg05_500K <- data.table::fread(input = "./obsDTg05_500K.csv", header = TRUE, na.strings = "NA_h2o")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  setkeyv(Odat_DT, cols = c("ID", "t"))
  head(Odat_DT)
  nrow(Odat_DT)

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

  # ----------------------------------------------------------------
  # IMPORT DATA
  # ----------------------------------------------------------------
  # options(stremr.verbose = TRUE)
  # set_all_stremr_options(fit.package = "glm", fit.algorithm = "glm")
  set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  # To inspect the input data.table:
  # OData$dat.sVar

  # ------------------------------------------------------------------
  # Fit propensity scores for Treatment, Censoring & Monitoring
  # ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))

  gform_CENS <- c("C ~ highA1c + t")
  # stratify_CENS <- list(C=c("t < 16", "t == 16"))
  # stratify_CENS <- list()

  gform_MONITOR <- "N ~ 1"
  # **** really want to define it like this ****
  # gform_TRT = c(list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t==0),
  #               list("TI[t] ~ CVD[t] + highA1c[t] + N[t-1]", t>0))

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  surv1 <- survNPMSM(wts.St.dlow, OData)

  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
  surv2 <- survNPMSM(wts.St.dhigh, OData)

  surv_by_trt <- list(surv1[['estimates']][["St.NPMSM"]], surv2[['estimates']][["St.NPMSM"]])
  names(surv_by_trt) <- c(surv1[['estimates']][["rule.name"]][1], surv2[['estimates']][["rule.name"]][1])
  t_idx <- surv1[['estimates']][['time']]

  f_plot_survest(surv_by_trt, t_idx)
  pl <- ggsurv(list(surv1[["estimates"]], surv2[["estimates"]]))


  # ------------------------------------------------------------------
  # Testing for bug with no surv at t=0 when no events at t=0 have occurred
  # ------------------------------------------------------------------
  wts.St.dlow_test <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
  wts.St.dlow_test[t==0, ("Y.tplus1") := 0]
  wts.St.dlow_test[t==1, ("Y.tplus1") := 0]
  wts.St.dlow_test2 <- wts.St.dlow_test[t > 0, ]
  surv1_test <- survNPMSM(wts.St.dlow_test2, OData)
  checkEquals(sum(surv1_test[["estimates"]][["time"]] == 0L)==1, TRUE)

  surv1_test2 <- survDirectIPW(wts.St.dlow_test, OData)
  checkEquals(sum(surv1_test[["estimates"]][["time"]] == 0L)==1, TRUE)



  # ------------------------------------------------------------------
  # Piping the workflow
  # ------------------------------------------------------------------
  require("magrittr")
  St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly") %>%
             survNPMSM(OData)  %$%
             estimates
  St.dlow
  #     t   sum_Y_IPAW sum_all_IPAW          ht   St.NPMSM      ht.KM     St.KM           rule.name
  # 1   0 3.978748e+01 1.298290e+03 0.030646065 0.9693539 0.03088578 0.9691142 gTI.dlowgPois3.yrly
  # 2   1 7.117673e+01 4.233691e+03 0.016811981 0.9530572 0.01142514 0.9580420 gTI.dlowgPois3.yrly
  # 3   2 1.228657e+02 1.291770e+04 0.009511416 0.9439923 0.02007299 0.9388112 gTI.dlowgPois3.yrly
  # 4   3 7.774297e+02 3.819757e+04 0.020352859 0.9247793 0.01303538 0.9265734 gTI.dlowgPois3.yrly
  # 5   4 2.923608e+03 1.185585e+05 0.024659625 0.9019746 0.01761006 0.9102564 gTI.dlowgPois3.yrly
  # 6   5 8.105229e+03 3.572692e+05 0.022686614 0.8815119 0.02240717 0.8898601 gTI.dlowgPois3.yrly
  # 7   6 1.731353e+04 1.145010e+06 0.015120852 0.8681826 0.01964637 0.8723776 gTI.dlowgPois3.yrly
  # 8   7 2.656852e+04 3.617676e+06 0.007344086 0.8618066 0.01603206 0.8583916 gTI.dlowgPois3.yrly
  # 9   8 1.326506e+04 1.218998e+07 0.001088193 0.8608688 0.01629328 0.8444056 gTI.dlowgPois3.yrly
  # 10  9 2.922616e+06 3.766485e+07 0.077595330 0.7940694 0.02484472 0.8234266 gTI.dlowgPois3.yrly
  # 11 10 1.692784e+06 6.843835e+07 0.024734435 0.7744286 0.02264685 0.8047786 gTI.dlowgPois3.yrly
  # 12 11 2.695594e+06 1.029670e+08 0.026179202 0.7541546 0.02172339 0.7872960 gTI.dlowgPois3.yrly
  # 13 12 1.438907e+06 1.388238e+08 0.010364986 0.7463378 0.01850481 0.7727273 gTI.dlowgPois3.yrly
  # 14 13 2.008508e+06 1.770313e+08 0.011345500 0.7378703 0.02488688 0.7534965 gTI.dlowgPois3.yrly
  # 15 14 1.822200e+06 2.105249e+08 0.008655511 0.7314836 0.01701469 0.7406760 gTI.dlowgPois3.yrly
  # 16 15 8.948983e+06 2.536594e+08 0.035279529 0.7056772 0.02202990 0.7243590 gTI.dlowgPois3.yrly
  # 17 16 0.000000e+00 0.000000e+00         NaN       NaN         NA        NA gTI.dlowgPois3.yrly

  St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow") %>%
             survNPMSM(OData)  %$%
             estimates
  St.dlow
#       t   sum_Y_IPAW sum_all_IPAW         ht   St.NPMSM      ht.KM     St.KM rule.name
# 1   0 2.497729e+01 9.771720e+02 0.02556079 0.9744392 0.03088578 0.9691142  gTI.dlow
# 2   1 2.514466e+01 2.450115e+03 0.01026265 0.9644389 0.01142514 0.9580420  gTI.dlow
# 3   2 1.314491e+02 6.420140e+03 0.02047450 0.9446925 0.02007299 0.9388112  gTI.dlow
# 4   3 2.036853e+02 1.664877e+04 0.01223425 0.9331349 0.01303538 0.9265734  gTI.dlow
# 5   4 7.413003e+02 4.423577e+04 0.01675794 0.9174975 0.01761006 0.9102564  gTI.dlow
# 6   5 2.443747e+03 1.174626e+05 0.02080447 0.8984094 0.02240717 0.8898601  gTI.dlow
# 7   6 5.864350e+03 3.118531e+05 0.01880485 0.8815150 0.01964637 0.8723776  gTI.dlow
# 8   7 1.279201e+04 8.393273e+05 0.01524078 0.8680800 0.01603206 0.8583916  gTI.dlow
# 9   8 3.544700e+04 2.295992e+06 0.01543864 0.8546780 0.01629328 0.8444056  gTI.dlow
# 10  9 1.535462e+05 6.333759e+06 0.02424251 0.8339585 0.02484472 0.8234266  gTI.dlow
# 11 10 3.833350e+05 1.732524e+07 0.02212582 0.8155065 0.02264685 0.8047786  gTI.dlow
# 12 11 1.035996e+06 4.803319e+07 0.02156833 0.7979173 0.02172339 0.7872960  gTI.dlow
# 13 12 2.353216e+06 1.346120e+08 0.01748148 0.7839686 0.01850481 0.7727273  gTI.dlow
# 14 13 8.985192e+06 3.836015e+08 0.02342324 0.7656055 0.02488688 0.7534965  gTI.dlow
# 15 14 1.730826e+07 1.089743e+09 0.01588289 0.7534455 0.01701469 0.7406760  gTI.dlow
# 16 15 2.800000e+07 1.271000e+09 0.02202990 0.7368471 0.02202990 0.7243590  gTI.dlow
# 17 16 0.000000e+00 0.000000e+00        NaN       NaN         NA        NA  gTI.dlow

  St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly") %>%
              survNPMSM(OData) %$%
              estimates
  St.dhigh

  #     t   sum_Y_IPAW sum_all_IPAW          ht   St.NPMSM       ht.KM     St.KM            rule.name
  # 1   0 2.289231e+02     38797.58 0.005900449 0.9940996 0.006918386 0.9930816 gTI.dhighgPois3.yrly
  # 2   1 3.421184e+03    167454.05 0.020430586 0.9737895 0.014900450 0.9782843 gTI.dhighgPois3.yrly
  # 3   2 1.975186e+04    658916.11 0.029976293 0.9445989 0.023289455 0.9555005 gTI.dhighgPois3.yrly
  # 4   3 9.358293e+04   2405269.41 0.038907461 0.9078470 0.026784229 0.9299082 gTI.dhighgPois3.yrly
  # 5   4 3.798779e+05   8145915.00 0.046634158 0.8655103 0.031530671 0.9005876 gTI.dhighgPois3.yrly
  # 6   5 1.394217e+06  27129914.41 0.051390409 0.8210314 0.033530572 0.8703904 gTI.dhighgPois3.yrly
  # 7   6 4.781450e+06  86031333.00 0.055578011 0.7754001 0.031303919 0.8431437 gTI.dhighgPois3.yrly
  # 8   7 1.485853e+07 224020681.36 0.066326577 0.7239704 0.026079870 0.8211546 gTI.dhighgPois3.yrly
  # 9   8 1.546769e+07 319653331.43 0.048388937 0.6889383 0.026720883 0.7992127 gTI.dhighgPois3.yrly
  # 10  9 1.497376e+07 374710197.48 0.039960904 0.6614077 0.023041475 0.7807976 gTI.dhighgPois3.yrly
  # 11 10 1.879632e+07 460587948.74 0.040809415 0.6344160 0.018837803 0.7660891 gTI.dhighgPois3.yrly
  # 12 11 1.135891e+07 542599644.97 0.020934244 0.6211350 0.014435696 0.7550301 gTI.dhighgPois3.yrly
  # 13 12 1.043877e+07 635755651.25 0.016419467 0.6109363 0.014784946 0.7438670 gTI.dhighgPois3.yrly
  # 14 13 1.597339e+07 728129607.15 0.021937557 0.5975338 0.015726496 0.7321686 gTI.dhighgPois3.yrly
  # 15 14 1.150659e+07 790411071.79 0.014557735 0.5888351 0.013555787 0.7222435 gTI.dhighgPois3.yrly
  # 16 15 1.322386e+07 880925228.28 0.015011334 0.5799959 0.020091646 0.7077324 gTI.dhighgPois3.yrly
  # 17 16 0.000000e+00         0.00         NaN       NaN          NA        NA gTI.dhighgPois3.yrly


  St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh") %>%
              survNPMSM(OData) %$%
              estimates
  St.dhigh
#     t   sum_Y_IPAW sum_all_IPAW          ht   St.NPMSM       ht.KM     St.KM rule.name
# 1   0 1.473714e+02 2.952989e+04 0.004990583 0.9950094 0.006918386 0.9930816 gTI.dhigh
# 2   1 1.334365e+03 9.618653e+04 0.013872679 0.9812060 0.014900450 0.9782843 gTI.dhigh
# 3   2 6.690061e+03 2.919185e+05 0.022917568 0.9587191 0.023289455 0.9555005 gTI.dhigh
# 4   3 2.257649e+04 8.315638e+05 0.027149436 0.9326904 0.026784229 0.9299082 gTI.dhigh
# 5   4 8.192048e+04 2.302602e+06 0.035577344 0.8995078 0.031530671 0.9005876 gTI.dhigh
# 6   5 2.464290e+05 6.129349e+06 0.040204751 0.8633433 0.033530572 0.8703904 gTI.dhigh
# 7   6 7.432322e+05 1.625016e+07 0.045736924 0.8238566 0.031303919 0.8431437 gTI.dhigh
# 8   7 1.940148e+06 4.260883e+07 0.045533949 0.7863432 0.026079870 0.8211546 gTI.dhigh
# 9   8 6.019081e+06 1.150543e+08 0.052315142 0.7452055 0.026720883 0.7992127 gTI.dhigh
# 10  9 1.471270e+07 3.044536e+08 0.048324925 0.7091935 0.023041475 0.7807976 gTI.dhigh
# 11 10 2.753646e+07 7.021941e+08 0.039214885 0.6813826 0.018837803 0.7660891 gTI.dhigh
# 12 11 2.490101e+07 1.332969e+09 0.018680861 0.6686538 0.014435696 0.7550301 gTI.dhigh
# 13 12 3.114277e+07 2.191312e+09 0.014211931 0.6591509 0.014784946 0.7438670 gTI.dhigh
# 14 13 4.006447e+07 2.773207e+09 0.014446980 0.6496282 0.015726496 0.7321686 gTI.dhigh
# 15 14 3.900000e+07 2.877000e+09 0.013555787 0.6408219 0.013555787 0.7222435 gTI.dhigh
# 16 15 5.700000e+07 2.837000e+09 0.020091646 0.6279468 0.020091646 0.7077324 gTI.dhigh
# 17 16 0.000000e+00 0.000000e+00         NaN       NaN          NA        NA gTI.dhigh


  # ------------------------------------------------------------------
  # Running IPW-adjusted MSM for the hazard
  # ------------------------------------------------------------------
  MSM.IPAW <- survMSM(OData,
                      wts_data = list(dlow = wts.St.dlow, dhigh = wts.St.dhigh),
                      t_breaks = c(1:8,12,16)-1,
                      est_name = "IPAW", getSEs = TRUE)

  MSM.IPAW[["estimates"]]
  attributes(MSM.IPAW[["estimates"]][[1]])
  pl <- ggsurv(MSM.IPAW[["estimates"]])
  # pl <- pl + ggplot2::theme(legend.position = "none")

  # names(MSM.IPAW)
  # nrow(t(MSM.IPAW$IC.Var.S.d$gTI.dlow$IC.S))
  # nrow(t(MSM.IPAW$IC.Var.S.d$gTI.dhigh$IC.S))
  # ncol(t(MSM.IPAW$IC.Var.S.d$gTI.dlow$IC.S))
  # ncol(t(MSM.IPAW$IC.Var.S.d$gTI.dhigh$IC.S))
  # # MSM.IPAW$St
  # names(MSM.IPAW$wts_data)
  # OData$nuniqueIDs

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP
  # ---------------------------------------------------------------------------------------------------------
  # t.surv <- c(4,5)
  t.surv <- c(9,10)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  # params = list(fit.package = "speedglm", fit.algorithm = "glm")
  # models = params,
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est3[["estimates"]][]
  #    est_name     t      risk      surv ALLsuccessTMLE nFailedUpdates   type rule.name
  # 1:    GCOMP    10 0.6600633 0.3399367          FALSE             11 pooled gTI.dhigh
  pl <- ggsurv(list(gcomp_est3[["estimates"]]))

  # models = params,
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  #    est_name     t       risk     surv ALLsuccessTMLE nFailedUpdates       type TMLE_Var  TMLE_SE rule.name
  # 1:     TMLE    10 0.04176398 0.958236           TRUE              0 stratified 690061.7 830.6995 gTI.dhigh
  pl <- ggsurv(list(tmle_est3[["estimates"]]))

  # pooling all observations (no stratification):
  # models = params,
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]
#    est_name     t       risk     surv ALLsuccessTMLE nFailedUpdates   type TMLE_Var  TMLE_SE rule.name
# 1:     TMLE    10 0.04172503 0.958275           TRUE              0 pooled 690061.7 830.6995 gTI.dhigh

  # ------------------------------------------------------------------------
  # RUN in PARALLEL seq-GCOMP & TMLE over t.surv (MUCH FASTER)
  # ------------------------------------------------------------------------
  # require("doParallel")
  # registerDoParallel(cores = 2)
  if (exists("setDTthreads")) data.table::setDTthreads(1)

  t.surv <- c(0,1,4)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  gcomp_est1 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "pooledGCOMP.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est1$estimates[]
# 1:    GCOMP     0 0.5018253 0.4981747          FALSE              1 pooled pooledGCOMP.dhigh
# 2:    GCOMP     1 0.6243751 0.3756249          FALSE              2 pooled pooledGCOMP.dhigh
# 3:    GCOMP     4 0.6597016 0.3402984          FALSE              5 pooled pooledGCOMP.dhigh

  gcomp_est2 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "pooledGCOMP.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est2$estimates[]
# 1:    GCOMP     0 0.5018090 0.4981910          FALSE              1 pooled pooledGCOMP.dlow
# 2:    GCOMP     1 0.6232914 0.3767086          FALSE              2 pooled pooledGCOMP.dlow
# 3:    GCOMP     4 0.6592231 0.3407769          FALSE              5 pooled pooledGCOMP.dlow

  # tmle_est_par1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "pool.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  tmle_est_par1 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "pooledTMLE.dhigh", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = FALSE)
  tmle_est_par1$estimates[]
# 1:     TMLE     0 0.005105404 0.9948946           TRUE              0 pooled 5.187611e-06 0.002277633 pooledTMLE.dhigh
# 2:     TMLE     1 0.018987340 0.9810127           TRUE              0 pooled 2.644945e-04 0.016263287 pooledTMLE.dhigh
# 3:     TMLE     4 0.041725030 0.9582750           TRUE              0 pooled 1.683729e+00 1.297585934 pooledTMLE.dhigh

  # tmle_est_par2 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "pool.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = TRUE)
  tmle_est_par2 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "pooledTMLE.dlow", Qforms = Qforms, stratifyQ_by_rule = FALSE, parallel = FALSE)
  tmle_est_par2$estimates[]
# 1:     TMLE     0 0.02517066 0.9748293           TRUE              0 pooled 3.918404e-06 0.001979496 pooledTMLE.dlow
# 2:     TMLE     1 0.03499178 0.9650082           TRUE              0 pooled 1.588953e-05 0.003986167 pooledTMLE.dlow
# 3:     TMLE     4 0.06577697 0.9342230           TRUE              0 pooled 8.724195e-03 0.093403400 pooledTMLE.dlow

  # tmle_est_par3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "strat.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE, parallel = TRUE)
  tmle_est_par3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", rule_name = "stratTMLE.dhigh", Qforms = Qforms, stratifyQ_by_rule = TRUE, parallel = FALSE)
  tmle_est_par3$estimates[]
# 1:     TMLE     0 0.005106787 0.9948932           TRUE              0 stratified 5.187316e-06 0.002277568 stratTMLE.dhigh
# 2:     TMLE     1 0.018998648 0.9810014           TRUE              0 stratified 2.644692e-04 0.016262509 stratTMLE.dhigh
# 3:     TMLE     4 0.041763985 0.9582360           TRUE              0 stratified 1.683719e+00 1.297581999 stratTMLE.dhigh

  # tmle_est_par4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "strat.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE, parallel = TRUE)
  tmle_est_par4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dlow", rule_name = "stratTMLE.dlow", Qforms = Qforms, stratifyQ_by_rule = TRUE, parallel = FALSE)
  tmle_est_par4$estimates[]
# 1:     TMLE     0 0.02504999 0.9749500           TRUE              0 stratified 3.916233e-06 0.001978948 stratTMLE.dlow
# 2:     TMLE     1 0.03482876 0.9651712           TRUE              0 stratified 1.588743e-05 0.003985904 stratTMLE.dlow
# 3:     TMLE     4 0.06549796 0.9345020           TRUE              0 stratified 8.724938e-03 0.093407376 stratTMLE.dlow

  pres <- ggsurv(list(gcomp_est1[["estimates"]], gcomp_est2[["estimates"]]))
  # pres

  pl <- ggsurv(list(tmle_est_par1[["estimates"]], tmle_est_par2[["estimates"]]), CI_line = TRUE)
  pl <- ggsurv(list(tmle_est_par1[["estimates"]], tmle_est_par2[["estimates"]]), CI_line = FALSE)
  pl <- ggsurv(list(tmle_est_par1[["estimates"]], tmle_est_par2[["estimates"]]), surv_col = c('red', 'black'))
  pl <- ggsurv(list(tmle_est_par1[["estimates"]], tmle_est_par2[["estimates"]]), surv_col = c('red', 'black'), CI_line = FALSE)

  # ------------------------------------------------------------------
  # Make a report:
  # ------------------------------------------------------------------
  # report.path <- "path/to/report/dir"
  # file.path = report.path,
  # test for opening file in local OS
  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, file.name = "sim.data.example.fup2", title = "Custom", author = "Insert Author Name",
                    openFile = TRUE)
                    # openFile = FALSE)

  # tmle_est_par1[["estimates"]]

  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  FUPtables = get_FUPtimes(MSM.IPAW$wts_data, IDnode = "ID", tnode = "t"),
                  openFile = TRUE,
                  # openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", x_legend = 9.5, y_legend = 0.99,
                  save_report_data = TRUE
                  )
  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  FUPtables = get_FUPtimes(MSM.IPAW$wts_data, IDnode = "ID", tnode = "t"),
                  openFile = TRUE,
                  # openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", x_legend = 9.5, y_legend = 0.99,
                  CI_line = FALSE
                  )

  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  FUPtables = get_FUPtimes(MSM.IPAW$wts_data, IDnode = "ID", tnode = "t"),
                  openFile = TRUE,
                  # openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", x_legend = 9.5, y_legend = 0.99,
                  use_ggplot = FALSE
                  )


  load(file = file.path(getOption('stremr.file.path'), "sim.data.example.fup") %+% ".Rd")
  names(report_results_list)

  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = report_results_list$NPMSM, MSM = report_results_list$MSM, GCOMP = report_results_list$GCOMP, TMLE = report_results_list$TMLE,
                  format = "html",
                  FUPtables = report_results_list$FUPtables,
                  # openFile = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = report_results_list$MSM.RDtables,
                  TMLE.RDtables = report_results_list$TMLE.RDtables,
                  WTtables = report_results_list$WTtables,
                  file.name = "sim.data.example.fup",
                  title = "Custom Report Title", author = "Insert Author Name", x_legend = 9.5, y_legend = 0.99
                  )

  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  AddFUPtables = TRUE,
                  # openFile = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name")

  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  AddFUPtables = TRUE,
                  # openFile = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", x_legend = "topright")


  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), plotKM = TRUE, MSM = MSM.IPAW, GCOMP = list(gcomp_est1, gcomp_est2), TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "html",
                  AddFUPtables = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", x_legend = 9.5, y_legend = 0.99)

  if (rmarkdown::pandoc_available(version = "1.12.3"))
    make_report_rmd(OData, NPMSM = list(surv1, surv2), MSM = MSM.IPAW, TMLE = list(tmle_est_par1, tmle_est_par2),
                  format = "pdf",
                  AddFUPtables = TRUE,
                  openFile = FALSE,
                  MSM.RDtables = get_MSM_RDs(MSM.IPAW, t.periods.RDs = c(12, 15), getSEs = TRUE),
                  TMLE.RDtables = get_TMLE_RDs(list(tmle_est_par1, tmle_est_par2), t.periods.RDs = c(1, 4)),
                  WTtables = get_wtsummary(MSM.IPAW$wts_data, cutoffs = c(0, 0.5, 1, 10, 20, 30, 40, 50, 100, 150), by.rule = TRUE),
                  file.name = "sim.data.example.fup", title = "Custom Report Title", author = "Insert Author Name", x_legend = 9.5, y_legend = 0.99)

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a stochastic intervention on MONITOR
  # *** make a separate function
  # **** +1 add evaluation of IPW & TMLE under static NDE
  # **** +2 add evaluation of IPW & TMLE under stochastic NDE intervention on monitoring
  # **** +3 add evaluation of IPW & TMLE under stochastic intervention on treatment
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]
}


# ---------------------------------------------------------------------------
# Test speedglm
# ---------------------------------------------------------------------------
test.speedglm.stochastic.TMLE.NDE.1Kdata <- function() {
  options(stremr.verbose = FALSE)
  options(width = 100)
  `%+%` <- function(a, b) paste0(a, b)
  require("data.table")
  data(OdatDT_10K)
  Odat_DT <- OdatDT_10K
  Odat_DT <- Odat_DT[ID %in% (1:5000), ]
  # define intervention on N as 0101010101...
  Odat_DT[, ("N.star.0101") := t%%2]
  setkeyv(Odat_DT, cols = c("ID", "t"))

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
  # ----------------------------------------------------------------
  # IMPORT DATA
  # ----------------------------------------------------------------
  set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm")
  OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
  # ------------------------------------------------------------------
  # Fit propensity scores for Treatment, Censoring & Monitoring
  # ------------------------------------------------------------------
  gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT <- list(
    TI=c("t == 0L",                                            # MODEL TI AT t=0
         "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
         "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
         "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
        ))
  gform_CENS <- c("C ~ highA1c + t")
  gform_MONITOR <- "N ~ 1"

  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
                          stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)

  # ---------------------------------------------------------------------------------------------------------
  # IPW-KM with stochastic intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "gPois3.yrly")
  surv1.stoch <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly")
  surv2.stoch <- survNPMSM(wts.St.dhigh, OData)

  # ---------------------------------------------------------------------------------------------------------
  # IPW-KM with static intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101")
  surv1.stat <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101")
  surv2.stat <- survNPMSM(wts.St.dhigh, OData)

  # ---------------------------------------------------------------------------------------------------------
  # IPW-KM with static intervention on MONITOR under NDE assumption
  # ---------------------------------------------------------------------------------------------------------
  wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv1.statNDE <- survNPMSM(wts.St.dlow, OData)
  wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", useonly_t_MONITOR = "N.star.0101 == 1")
  surv2.statNDE <- survNPMSM(wts.St.dhigh, OData)

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a stochastic intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "gPois3.yrly", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a static intervention on MONITOR
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]

  # ---------------------------------------------------------------------------------------------------------
  # TMLE / GCOMP with a static intervention on MONITOR under NDE assumption
  # ---------------------------------------------------------------------------------------------------------
  t.surv <- c(4,5)
  Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
  params = list(fit.package = "speedglm", fit.algorithm = "glm")
  gcomp_est3 <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                            useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  gcomp_est3$estimates[]
  # stratified modeling by rule followers only:
  tmle_est3 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = TRUE)
  tmle_est3$estimates[]
  # pooling all observations (no stratification):
  tmle_est4 <- fitTMLE(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", intervened_MONITOR = "N.star.0101",
                        useonly_t_MONITOR = "N.star.0101 == 1", Qforms = Qforms, stratifyQ_by_rule = FALSE)
  tmle_est4$estimates[]
}
