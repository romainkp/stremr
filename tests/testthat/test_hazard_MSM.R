  require("data.table")
  require("stremr")
  options(stremr.verbose = FALSE)
  options(gridisl.verbose = FALSE)


  # ------------------------------------------------------------------------------------------------------
  # (IA) Data from the simulation study
  # ------------------------------------------------------------------------------------------------------
  # Nsize <- 1000
  # OdataNoCENS <- simulateDATA.fromDAG(Nsize = Nsize, rndseed = 124356)
  data(OdataNoCENS)
  OdataDT <- as.data.table(OdataNoCENS, key=c(ID, t))

  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataDT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]

  # Define intervention (always treated):
  OdataDT[, ("TI.set1") := 1L]
  OdataDT[, ("TI.set0") := 0L]

  gform_CENS <- "C ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  gform_MONITOR <- "N ~ 1"
  stratify_CENS <- list(C=c("t < 16", "t == 16"))
  OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "N.tminus1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  OData <- fitPropensity(OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, gform_MONITOR = gform_MONITOR, stratify_CENS = stratify_CENS)

  # ----------------------------------------------------------------------
  # IPW-MSM for hazard
  # ----------------------------------------------------------------------
  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
  ## this is your \bar{a} identifier of the static regimen 1, i.e., f(\bar{a}) create your own variable
  wts.DT.1[, "theta" := 1]
  wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
  ## this is your \bar{a} identifier of the static regimen 2, i.e., f(\bar{a}) create your own variable
  wts.DT.0[, "theta" := 0]

  ## old MSM that allows pooling over time intervals, theta defined above wont be used
  survMSM_res_old <- survMSM(list(wts.DT.1, wts.DT.0), OData, glm_package = "speedglm", return_wts = TRUE)

  survMSM_res_old[["TI1"]]

  ## new MSM that allows arbitrary pooling over time intervals, theta defined above can be used in the formula:
  survMSM_res <- survMSM2(list(wts.DT.1, wts.DT.0), form = "Y.tplus1 ~ t + theta + t:theta", OData)
  survMSM_res <- survMSM2(list(wts.DT.1, wts.DT.0), form = "Y.tplus1 ~ -1 + as.factor(t):as.factor(rule.name)", OData)
  survMSM_res[[1]]
  as.formula("Y.tplus1 ~ -1 + as.factor(t):as.factor(rule.name)")
  # survMSM_res$St

  str(terms(form))

  attributes(terms(form))$variables[[2]]
  attributes(terms(form))$variables[[3]]
  attributes(terms(form))$variables[[4]]

  labels(terms(form))

  newdata <- data.frame(t = 1:10, theta = 7.0)
  new_design_mat <- model.matrix(as.formula(form2), data = newdata)
  class(new_design_mat)


