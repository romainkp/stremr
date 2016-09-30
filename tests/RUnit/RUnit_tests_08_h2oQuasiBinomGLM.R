notest.h2oQuasiBinomGLM.Ensemble <- function() {
    reqh2o <- requireNamespace("h2o", quietly = TRUE)
    if (reqh2o) {
        options(stremr.verbose = FALSE)
        `%+%` <- function(a, b) paste0(a, b)
        require("h2o")
        h2o::h2o.init(nthreads = 1)
        # h2o::h2o.init(nthreads = -1)
        require('h2oEnsemble')
        require("data.table")
        # require("stremr")
        data(OdatDT_10K)
        Odat_DT <- OdatDT_10K

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


        OData <- importData(Odat_DT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "lastNat1.factor"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = outcome)
        # OData$define_CVfolds(nfolds = 10, seed = 23)

        # ---------------------------------------------------------------------------------------------------------
        # VALIDATING QUASIBINOMIAL (cont Y) LOGISTIC REG in H2O glm WITH TMLE
        # ---------------------------------------------------------------------------------------------------------
        t.surv <- c(10)
        Qforms <- rep.int("Q.kplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))
        # h2o::h2o.init(nthreads = -1, startH2O = FALSE)
        # h2o::h2o.shutdown()
        # options(stremr.verbose = TRUE)

        # speedglm:
        params = list(fit.package = "speedglm", fit.algorithm = "glm")
        gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
        gcomp_est$estimates[]
        #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates rule.name
        # 1:    GCOMP 10 0.2723941 0.7276059          FALSE             11 gTI.dhigh

        # H2O glm w/ L_BFGS:
        params = list(fit.package = "h2o", fit.algorithm = "glm", solver = "L_BFGS")
        gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
        gcomp_est$estimates[]
        #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates rule.name
        # 1:    GCOMP 10 0.2723858 0.7276142          FALSE             11 gTI.dhigh

        # H2O glm w/ IRLSM:
        params = list(fit.package = "h2o", fit.algorithm = "glm", solver = "IRLSM")
        gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = params, stratifyQ_by_rule = FALSE)
        gcomp_est$estimates[]
        #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates rule.name
        # 1:    GCOMP 10 0.2723941 0.7276059          FALSE             11 gTI.dhigh

        # H2O SL:
        h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
        h2o.glm.2 <- function(..., x = "highA1c", alpha = 0.0) h2o.glm.wrapper(..., x = x, alpha = alpha)
        glm_hyper_params <- list(search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
                                 alpha = c(0,1,seq(0.1,0.9,0.1)), lambda = c(0,1e-7,1e-5,1e-3,1e-1))
        SLparams <- list( fit.package = "h2o",
                         fit.algorithm = "SuperLearner",
                         grid.algorithm = c("glm"),
                         learner = c("h2o.glm.2"),
                         nfolds = 5,
                         seed = 23,
                         glm = glm_hyper_params)

        gcomp_est <- fitSeqGcomp(OData, t_periods = t.surv, intervened_TRT = "gTI.dhigh", Qforms = Qforms, params_Q = SLparams, stratifyQ_by_rule = FALSE)
        gcomp_est$estimates[]
        #    est_name  t      risk      surv ALLsuccessTMLE nFailedUpdates   type rule.name
        # 1:    GCOMP 10 0.2495067 0.7504933          FALSE             11 pooled gTI.dhigh

        h2o::h2o.shutdown(prompt = FALSE)
    }
}
