
#' @export
fitSeqGcomp <- function(OData,
                        t = OData$max.t,
                        Qforms,
                        stratify = FALSE,
                        params = list(),
                        verbose = getOption("stremr.verbose")) {


  gvars$verbose <- verbose
  nodes <- OData$nodes
  new.factor.names <- OData$new.factor.names
  assert_that(is.list(params))
  # ------------------------------------------------------------------------------------------------
  # Process the input formulas and stratification settings;
  # Define regression classes for Q.Y and put them in a single list of regressions.
  # ------------------------------------------------------------------------------------------------
  Qforms.default <- "Q.kplus1 ~ Lnodes + Anodes + Cnodes + Nnodes"
  stratify_Q <- list(nodes[['tnode']] %+% " == " %+% (OData$min.t:t))
  names(stratify_Q) <- "Q.kplus1"
  # names(stratify_Q[[1]]) <- rep.int("Q.kplus1", length(stratify_Q[[1]]))
  # ------------------------------------------------------------------------------------------------
  # TMLE:
  # Define a new column Q.kplus1.regimen for each regimen/trt of intersect
  # That column keeps the tabs on the running Q fit (sequantial g-comp)
  # ------------------------------------------------------------------------------------------------

  Q.sVars <- process_regforms(regforms = Qforms, default.reg = Qforms.default, stratify.EXPRS = stratify_Q, model_contrl = params,
                              OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = FALSE, outvar.class = "binary")

  Q.sVars$regs
  # g_CAN_regs_list[["Q.k"]] <- Q.sVars$regs
  # g_CAN_regs_list <- vector(mode = "list", length = 3)
  # names(g_CAN_regs_list) <- c("gC", "gA", "gN")
  # class(g_CAN_regs_list) <- c(class(g_CAN_regs_list), "ListOfRegressionForms")


  browser()

  # ------------------------------------------------------------------------------------------
  # DEFINE a single regression class
  # Perform S3 method dispatch on ALL_g_regs, which will determine the nested tree of SummaryModel objects
  # Perform fit and prediction
  # ------------------------------------------------------------------------------------------
  ALL_g_regs <- SGCompRegClass$new(RegressionForms = Q.sVars$regs)
  ALL_g_regs$S3class <- "generic"
  modelfits.g0 <- newsummarymodel(reg = ALL_g_regs, DataStorageClass.g0 = OData)
  modelfits.g0$fit(data = OData, predict = TRUE)
  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0
  h_gN <- modelfits.g0$predictAeqa(n = OData$nobs)

  # ------------------------------------------------------------------------------------------
  # Observed likelihood of (A,C,N) at each t, based on fitted object models in object modelfits.g0
  # ------------------------------------------------------------------------------------------
  # get back g_CAN_regs_list:
  OData$modelfits.g0 <- modelfits.g0

  ALL_g_regs <- modelfits.g0$reg
  g_CAN_regs_list <- ALL_g_regs$RegressionForms

  OData$modelfit.gC <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gC")]]
  OData$modelfit.gA <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gA")]]
  OData$modelfit.gN <- modelfits.g0$getPsAsW.models()[[which(names(g_CAN_regs_list) %in% "gN")]]

  g0.A <- OData$modelfit.gA$getcumprodAeqa()
  g0.C <- OData$modelfit.gC$getcumprodAeqa()
  g0.N <- OData$modelfit.gN$getcumprodAeqa()

  OData$dat.sVar[, c("g0.A", "g0.C", "g0.N", "g0.CAN") := list(g0.A, g0.C, g0.N, g0.A*g0.C*g0.N)]
  # newdat <- OData$dat.sVar[, list("g0.A" = g0.A, "g0.C" = g0.C, "g0.N" = g0.N, "g0.CAN" = g0.A*g0.C*g0.N)]

  return(OData)
}