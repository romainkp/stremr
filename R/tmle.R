# newsummarymodel.SGCompRegClass <- function(reg, DataStorageClass.g0, ...) GenericModel$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
# SGCompRegClass <- R6Class("SGCompRegClass",
#   inherit = RegressionClass,
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     # outvar = character(),          # the outcome variable name (Ynode)
#     # outvar.class = character(),
#     Qforms = character(),
#     # n_regs = integer(),
#     # n_timepts = integer(),       # number of time points
#     # nodes = list(),
#     # predvars = list(),             # list of predictors for each regression from Qforms
#     # not used for now:
#     Anodes = list(),               # list of treatment var names, by timepoint
#     Cnodes = list(),               # list of censoring var names, by timepoint
#     Lnodes = list(),               # list of time-varying confounder names, by timepoint
#     Wnodes = character(),          # character vector of baseline covariate names
#     subset = NULL,                 # subset expression (later evaluated to logical vector in the envir of the data)
#     ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
#     pool = logical()
#   )
# )

# SingleQLearnRegClass <- R6Class("SingleQLearnRegClass",
#   inherit = SingleRegressionFormClass,
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     reg_idx = integer(),
#     subset_censored = logical(),
#     stratify = TRUE,
#     pool_regimes = FALSE
#   )
# )

RegressionClassQlearn <- R6Class("RegressionClassQlearn",
  inherit = RegressionClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    reg_idx = integer(),
    subset_censored = logical(),
    stratify = TRUE,
    pool_regimes = FALSE,
    initialize = function(reg_idx, subset_censored, stratify, pool_regimes, ...) {
      if (!missing(reg_idx)) self$reg_idx <- reg_idx
      if (!missing(subset_censored)) self$subset_censored <- subset_censored
      if (!missing(stratify)) self$stratify <- stratify
      if (!missing(pool_regimes)) self$pool_regimes <- pool_regimes
      super$initialize(...)
    }
  ),
  active = list(
    get.reg = function() {
      list(reg_idx = self$reg_idx,
           outvar = self$outvar,
           predvars = self$predvars,
           outvar.class = self$outvar.class,
           subset_vars = self$subset_vars,
           subset_exprs = self$subset_exprs,
           subset_censored = self$subset_censored,
           stratify = self$stratify,
           pool_regimes = self$pool_regimes,
           model_contrl = self$model_contrl,
           censoring = self$censoring
           )
    }
  )
)

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
  # **** The stratification by follow-up has to be based only on 't' values that were observed in the data****
  # ------------------------------------------------------------------------------------------------
  Qperiods <- rev(OData$min.t:t)
  Qreg_idx <- rev(seq_along(Qperiods))
  stratify_Q <- as.list(nodes[['tnode']] %+% " == " %+% (Qperiods))
  names(stratify_Q) <- rep.int("Q.kplus1", length(stratify_Q))
  # ------------------------------------------------------------------------------------------------
  # Process the input formulas and stratification settings;
  # Define regression classes for Q.Y and put them in a single list of regressions.
  # ------------------------------------------------------------------------------------------------
  Qforms.default <- rep.int("Q.kplus1 ~ Lnodes + Anodes + Cnodes + Nnodes", length(Qperiods))
  # ------------------------------------------------------------------------------------------------
  # TMLE:
  # Initiate the Q.kplus1 - need to do this for each regimen
  # That column keeps the tabs on the running Q fit (SEQ G-COMP)
  OData$dat.sVar[, "Q.kplus1" := as.numeric(get(OData$nodes$Ynode))]
  OData$def.types.sVar()
  # ------------------------------------------------------------------------------------------------

  if (missing(Qforms)) Qforms <- Qforms.default

  Q_regs_list <- vector(mode = "list", length = length(stratify_Q))
  names(Q_regs_list) <- unlist(stratify_Q)
  class(Q_regs_list) <- c(class(Q_regs_list), "ListOfRegressionForms")

  for (i in seq_along(Q_regs_list)) {
    regform <- process_regform(as.formula(Qforms[[i]]), sVar.map = nodes, factor.map = new.factor.names)
    # old version:
    # Q.sVars <- process_regforms(regforms = Qforms, default.reg = Qforms.default, stratify.EXPRS = stratify_Q, model_contrl = params,
    #                             OData = OData, sVar.map = nodes, factor.map = new.factor.names, censoring = FALSE, outvar.class = "Qlearn")
    reg <- RegressionClassQlearn$new(reg_idx = Qreg_idx[i], outvar = "Q.kplus1", predvars = regform$predvars, outvar.class = list("Qlearn"),
                                     subset_vars = list("Q.kplus1"), subset_exprs = stratify_Q[i], model_contrl = params,
                                     censoring = FALSE)
    Q_regs_list[[i]] <- reg
  }
  Qlearn.fit <- GenericModel$new(reg = Q_regs_list, DataStorageClass.g0 = OData)

  browser()

  # ------------------------------------------------------------------------------------------
  # Note: since new regimen results in new Q.kplus1 and hence new outcome -> new regression, there is no point in doing all regimens at once.
  # Need to define A^* at the begining, as a column (one for each regimen)
  # ------------------------------------------------------------------------------------------
  # Need to implement:
  # ------------------------------------------------------------------------------------------
  # *** Need to add to current subsets expression for conditioning on uncensored observation (or censored) at each t
  # Q_regs_list_2 <- stratify_by_uncensored(Q_regs_list)
  # OData$dat.sVar[t == 0L, ]
  # ------------------------------------------------------------------------------------------
  # 0. CONSIDER THE PROS AND CONS OF KEEPING THE OLD CLASS SYSTEM (RegressionClass -> SingleRegressionClass) VS. STARTING FROM SCRATCH
  #    Another alternative is to just simply add to model controls whatever new information we need? And then keep everything indact
  #    Would require new model controls for each Q regression (with idx/subsets/etc)
  # 1. At each t iteration:
  #   Fit: Remove all censored at t when fitting (add new expr to subset defn or pass second subset expr?).
  #        Outcome is always Q.kplus1 (either by regimen or not).
  #   Prediction: For each Q.m (Q fit), rename columns to replace current A[t] with A^*[t], (possibly for each regimen), put A[t] back afterwards.
  #   Prediction: Add observation which were censored at t (C[t]==1) to the prediction set, using the second subset
  #   Prediction: Save ProbA1 from PredictP1() for all obs used in prediction in rows of column Q.kplus1 (either at row Q.kplus1[t] or at Q.kplus1[t-1])
  # ------------------------------------------------------------------------------------------
  # 2. At next iteration t-1:
  #   Fit: If previous prediction was saved at Q.kplus1[t-1], then all outcomes are already set to needed values
  #        If not, the outcomes at Q.kplus1[t-1] need to be updated with values from Q.kplus1[t] for all IDs who were used in prediction at t.
  #
  # ------------------------------------------------------------------------------------------
  # WIll run all regressions (one for each subsets defined above):
  Qlearn.fit$fit(data = OData, predict = TRUE)

  OData$dat.sVar[1:50,]

  Qlearn.fit
  # get the joint likelihood at each t for all 3 variables at once (P(C=c|...)P(A=a|...)P(N=n|...)).
  # NOTE: Separate predicted probabilities (e.g., P(A=a|...)) are also stored in individual child classes.
  # They are accessed later from modelfits.g0
  h_gN <- Qlearn.fit$predictAeqa(n = OData$nobs)
  class(Qlearn.fit$getPsAsW.models()[[1]]$getPsAsW.models()[[1]])



  # ------------------------------------------------------------------------------------------------
  # Doing the same as above but directly via SingleRegressionFormClass/SingleQLearnRegClass
  # ------------------------------------------------------------------------------------------------
  if (missing(Qforms)) {
    Qforms <- Qforms.default
  }
  stratify_Q <- list(Q.kplus1 = nodes[['tnode']] %+% " == " %+% (Qperiods))
  res <- process_regform(as.formula(Qforms[[1]]), sVar.map = nodes, factor.map = new.factor.names)
  reg <- SingleRegressionFormClass$new(outvar = "Q.kplus1", predvars = res$predvars, outvar.class = list("Qlearn"),
                                          subset_vars = list("Q.kplus1"), subset_exprs = stratify_Q, model_contrl = params,
                                          censoring = FALSE)

  ALL_Q_regs <- RegressionClass$new(RegressionForms = Q_regs_list)
  ALL_Q_regs$S3class <- "generic"
  class(ALL_Q_regs)
  # Qlearn.fit <- newsummarymodel(reg = ALL_Q_regs, DataStorageClass.g0 = OData)

  return(OData)
}