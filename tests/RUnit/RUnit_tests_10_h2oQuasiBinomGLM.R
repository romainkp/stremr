# notest.h2oQuasiBinomGLM.Ensemble <- function() {
#     reqh2o <- requireNamespace("h2o", quietly = TRUE)
#     if (reqh2o) {
#         options(stremr.verbose = TRUE)
#         `%+%` <- function(a, b) paste0(a, b)
#         require("h2o")
#         h2o::h2o.init(nthreads = 1)
#         # h2o::h2o.shutdown(prompt = FALSE)
#         # h2o::h2o.init(nthreads = -1)
#         # require('h2oEnsemble')
#         require("data.table")
#         # require("stremr")
#         data(OdatDT_10K)
#         Odat_DT <- OdatDT_10K

#         # ---------------------------------------------------------------------------
#         # Define some summaries (lags C[t-1], A[t-1], N[t-1])
#         # ---------------------------------------------------------------------------
#         ID <- "ID"; t <- "t"; TRT <- "TI"; I <- "highA1c"; outcome <- "Y.tplus1";
#         lagnodes <- c("C", "TI", "N")
#         newVarnames <- lagnodes %+% ".tminus1"
#         Odat_DT[, (newVarnames) := shift(.SD, n=1L, fill=0L, type="lag"), by=ID, .SDcols=(lagnodes)]
#         # indicator that the person has never been on treatment up to current t
#         Odat_DT[, ("barTIm1eq0") := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
#         Odat_DT[, ("lastNat1.factor") := as.factor(lastNat1)]

#     # ------------------------------------------------------------------
#     # Propensity score models for Treatment, Censoring & Monitoring
#     # ------------------------------------------------------------------
#     gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"
#     stratify_TRT <- list(
#       TI=c("t == 0L",                                            # MODEL TI AT t=0
#            "(t > 0L) & (N.tminus1 == 1L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN MONITORED
#            "(t > 0L) & (N.tminus1 == 0L) & (barTIm1eq0 == 1L)",  # MODEL TRT INITATION WHEN NOT MONITORED
#            "(t > 0L) & (barTIm1eq0 == 0L)"                       # MODEL TRT CONTINUATION (BOTH MONITORED AND NOT MONITORED)
#           ))

#     gform_CENS <- c("C ~ highA1c + t")
#     gform_MONITOR <- "N ~ 1"


#     OData <- stremr::importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
#     # OData$define_CVfolds(nfolds = 10, seed = 23)
#     OData$fast.load.to.H2O(saveH2O = TRUE)
#     h2o.no_progress()
#     # h2o.show_progress()

#     # params = list(fit.package = "h2o", fit.algorithm = "glm", solver = "L_BFGS", family = "quasibinomial")
#     # set_all_stremr_options(fit.package = "h2o", fit.algorithm = "glm") # TO DO: add family to glob options , family = "binomial"
#     # set_all_stremr_options(fit.package = "speedglm", fit.algorithm = "glm") # TO DO: add family to glob options , family = "binomial"

#     stremrOptions("fit.package", "h2o")
#     stremrOptions("fit.algorithm", "glm")
#     model <- "h2o.glm"

#     params_CENS = list(ntrees = 5000, learn_rate = 0.01, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)
#     params_TRT = list(ntrees = 5000, learn_rate = 0.01, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)
#     params_MONITOR = list(ntrees = 5000, learn_rate = 0.01, sample_rate = 0.8, col_sample_rate = 0.8, balance_classes = TRUE)

#     OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT,
#                             stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR,
#                             params_CENS = params_CENS, params_TRT = params_TRT, params_MONITOR = params_MONITOR)

#     wts.St.dlow <- getIPWeights(OData, intervened_TRT = "gTI.dlow")
#     surv1 <- survNPMSM(wts.St.dlow, OData)
#     wts.St.dhigh <- getIPWeights(OData, intervened_TRT = "gTI.dhigh")
#     surv2 <- survNPMSM(wts.St.dhigh, OData)


#   #       h2oFrame <- as.h2o(Odat_DT)
#   #       h2oFrame[, "Qkplus1"] <-  h2o.asnumeric(h2oFrame[["Y.tplus1"]]) # set the initial values of Q (the observed outcome node)
#   #       h2oFrame[, "prev_Qkplus1"] <- h2oFrame[, "Qkplus1"]
#   #       set.seed(435)
#   #       h2oFrame <- h2o::h2o.createFrame()
#   #       set.seed(435)
#   #       h2oFrame2 <- h2o::h2o.createFrame(randomize = FALSE)

#   #       subset_idx <- sort(sample(1:nrow(h2oFrame), as.integer(nrow(h2oFrame)/3)))

#   #       h2oFrame[subset_idx,]

#   #       as.data.table(h2oFrame[h2oFrame[["prev_Qkplus1"]] != h2oFrame[["Qkplus1"]], ])

#   #       prev_Qkplus1 <- h2oFrame[subset_idx, "Qkplus1"]
#   #       h2oFrame[subset_idx, "prev_Qkplus1"] <- prev_Qkplus1
#   #       as.data.table(h2oFrame[h2oFrame[["prev_Qkplus1"]] != h2oFrame[["Qkplus1"]], ])

#   # OData$H2Oframe[, "Qkplus1"] <-  h2o.asnumeric(OData$H2Oframe[[OData$nodes$Ynode]]) # set the initial values of Q (the observed outcome node)
#   # OData$H2Oframe[, "prev_Qkplus1"] <- OData$H2Oframe[, "Qkplus1"]

#         t.surv <- c(10)
#         Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

#         # H2O glm w/ L_BFGS:
#         params = list(fit.package = "h2o", fit.algorithm = "glm", solver = "L_BFGS", family = "quasibinomial")
#         gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params,
#                                  stratifyQ_by_rule = FALSE, TMLE = TRUE)
#         gcomp_est$estimates[]

#         ## 1. When swapping the treatment nodes to use correct counterfactual exposure for Q
#         ## 2. Making sure that row slicing is accomplished as intended:
#         ## SEQUENTIAL GCOMP:
#         # [1] "Surv est 1: 0.727554692299103"
#         # [1] "Surv est 2: 0.727587295034153"
#         ## TMLE:
#         # [1] "Surv est 1: 0.699256088690996"
#         # [1] "Surv est 2: 0.699290661377025"


#         ## When ignoring possibly incorrect subset (row slice) assigning:
#         # [1] "Surv est 1: 0.742087351492688"
#         # [1] "Surv est 2: 0.742088692289334"
#         ## When correcting and making sure the slice assignment is correct:
#         # [1] "Surv est 1: 0.742047150582591"
#         # [1] "Surv est 2: 0.742043225574081"

#         # PREVIOUS VALIDATION RUN:
#         #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates rule.name
#         # 1:    GCOMP 10 0.2723858 0.7276142          FALSE             11 gTI.dhigh



#         # ---------------------------------------------------------------------------------------------------------
#         # VALIDATING QUASIBINOMIAL (cont Y) LOGISTIC REG in H2O glm WITH TMLE
#         # ---------------------------------------------------------------------------------------------------------
#         # h2o::h2o.init(nthreads = -1, startH2O = FALSE)
#         # h2o::h2o.shutdown()
#         # options(stremr.verbose = TRUE)

#         # speedglm:
#         params = list(fit.package = "speedglm", fit.algorithm = "glm")
#         gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)
#         gcomp_est$estimates[]
#         #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates rule.name
#         # 1:    GCOMP 10 0.2723941 0.7276059          FALSE             11 gTI.dhigh

#         # H2O glm w/ IRLSM:
#         params = list(fit.package = "h2o", fit.algorithm = "glm", solver = "IRLSM", family = "quasibinomial")
#         gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)
#         gcomp_est$estimates[]
#         #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates rule.name
#         # 1:    GCOMP 10 0.2723941 0.7276059          FALSE             11 gTI.dhigh

#         # H2O SL:
#         h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
#         h2o.glm.2 <- function(..., x = "highA1c", alpha = 0.0) h2o.glm.wrapper(..., x = x, alpha = alpha)
#         glm_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
#                                  alpha = c(0,1,seq(0.1,0.9,0.1)), lambda = c(0,1e-7,1e-5,1e-3,1e-1))
#         SLparams <- list( fit.package = "h2o",
#                          fit.algorithm = "SuperLearner",
#                          grid.algorithm = c("glm"),
#                          learner = c("h2o.glm.2"),
#                          nfolds = 5,
#                          seed = 23,
#                          glm = glm_hyper_params)

#         gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, models = SLparams, stratifyQ_by_rule = FALSE)
#         gcomp_est$estimates[]
#         #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates   type rule.name
#         # 1:    GCOMP 10 0.2495067 0.7504933          FALSE             11 pooled gTI.dhigh

#         h2o::h2o.shutdown(prompt = FALSE)
#     }
# }
