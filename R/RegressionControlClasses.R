## ---------------------------------------------------------------------------------
## SPECIFYING regressions
## ---------------------------------------------------------------------------------
get_vars_fromlist <- function(varname, sVar.map) {
  if (varname %in% names(sVar.map)) {
    as.vector(sVar.map[[varname]])
  } else {
    varname
  }
}
# Parse the formulas for summary measure names and create a map to actual covariate names in sA & sW
process_regform <- function(regform, sVar.map = NULL, factor.map = NULL) {
  # regform1 <- as.formula("N ~ 1")
  # Getting predictors (sW names):
  regformterms <- terms(regform)
  sW.names <- attributes(regformterms)$term.labels
  sW.names.alt <- colnames(attributes(regformterms)$factors)
  assert_that(all(sW.names == sW.names.alt))
  # Getting OUTCOMEs (sA names):
  (out.var <- deparse(attributes(regformterms)$variables[[2]])) # LHS character string

  out.vars.form <- as.formula(". ~ " %+% out.var)
  out.vars.terms <- terms(out.vars.form)
  sA.names <- attributes(out.vars.terms)$term.labels

  outvars <- unlist(lapply(sA.names, get_vars_fromlist, sVar.map))
  predvars <- unlist(lapply(sW.names, get_vars_fromlist, sVar.map))
  # in case some factors were also involved (these will then be replaced only on a second iteration)
  predvars <- unlist(lapply(predvars, get_vars_fromlist, factor.map))
  return(list(outvars = outvars, predvars = predvars))
}

# Loop through a list of SingleRegressionFormClass objects and their outvars as if it was one long list of outvars and create the subsetting expressions
# This uses S3 method dispatch on object ListOfRegressionForms
stratify_by_uncensored <- function(regs) {
  for (Var_indx in seq_along(get_outvars(regs)[-1])) {
    curr_outvar <- get_outvars(regs)[Var_indx+1]
    curr_exprs <- get_subset_exprs(regs)[[Var_indx+1]]
    prev_outvars <- as.vector(unique(get_outvars(regs)[1:Var_indx]))
    prev_outvars_to_condition <- prev_outvars[!(prev_outvars %in% curr_outvar)]
    if (length(prev_outvars_to_condition) > 0) {
      strat.C <- paste0(prev_outvars_to_condition %+% " == " %+% gvars$noCENScat, collapse=" & ")
      if (!is.null(curr_exprs)) {
        reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_exprs = paste0(curr_exprs, " & ", strat.C))
        # reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_exprs = stringr::str_c(curr_exprs, " & ", strat.C))
      } else {
        reg.obj <- set_subset_exprs(regs, idx = Var_indx + 1, subset_exprs = strat.C)
      }
    }
  }
  return(regs)
}

# Create subsetting expressions for a node (Anode, Cnode or Nnode)
# Named list with character expressions for subsetting.
# Each list item corresponds to one outcome in SingleRegressionFormClass
create_subset_expr <- function(outvars, stratify.EXPRS) {
  if (is.null(stratify.EXPRS)) {
    return(NULL)
  }
  Node_subset_expr <- vector(mode="list", length=length(outvars))
  names(Node_subset_expr) <- outvars
  assert_that(is.list(stratify.EXPRS))
  if (!all(outvars %in% names(stratify.EXPRS))) {
    stop("Could not locate the appropriate regression variable(s) within the supplied stratification list stratify_CENS, stratify_TRT or stratify_MONITOR." %+% "\n" %+%
          "The regression outcome variable(s) specified in gform_CENS, gform_TRT or gform_MONITOR were: ( '" %+% paste0(outvars, collapse=",") %+% "' )" %+% "\n" %+%
          "However, the item names in the matching stratification list were: ( '" %+% paste0(names(stratify.EXPRS), collapse=",") %+% "' )"
          )
  }
  for (idx in seq_along(Node_subset_expr))
    if (!is.null(stratify.EXPRS[[outvars[idx]]]))
      Node_subset_expr[[idx]] <- stratify.EXPRS[[outvars[idx]]]
  return(Node_subset_expr)
}

# When several reg forms are specified (multivariate Anodes), process outvars into one vector and process predvars in a named list of vectors
# stratify.EXPRS - Must be a named list. One item (characeter vectors) per one outcome in regforms.
process_regforms <- function(regforms, default.reg, stratify.EXPRS = NULL, model_contrl = NULL, OData, sVar.map = NULL, factor.map = NULL, censoring = FALSE, outvar.class = NULL) {
  using.default <- FALSE
  if (missing(regforms)) {
    using.default <- TRUE
    regforms <- default.reg
  }
  if (!is.null(stratify.EXPRS)) assert_that(is.list(stratify.EXPRS))
  regs <- vector(mode="list", length=length(regforms))
  for (idx in seq_along(regforms)) {
    res <- process_regform(as.formula(regforms[[idx]]), sVar.map = sVar.map, factor.map = factor.map)
    if (using.default && gvars$verbose && !is.null(res$outvars))
      message("Using the default regression formula: " %+% paste0(res$outvars, collapse="+") %+% " ~ " %+% paste0(res$predvars, collapse="+"))

      if (!is.null(outvar.class)) {
        outvar.class <- as.list(rep.int(outvar.class, length(res$outvars)))
        names(outvar.class) <- res$outvars
      } else {
        # outvar.class <- rep.int(list("univariate"), length(res$outvars)) ## waz before now we automatically detect the outcome variable type
        outvar.class <- OData$type.sVar(res$outvars) 
        names(outvar.class) <- res$outvars
      }
      subset_exprs <- create_subset_expr(outvars = res$outvars, stratify.EXPRS = stratify.EXPRS)

      regobj <- RegressionClass$new(outvar = res$outvars, predvars = res$predvars, outvar.class = outvar.class,
                                    subset_vars = NULL, subset_exprs = subset_exprs, model_contrl = model_contrl,
                                    censoring = censoring)
      regs[[idx]] <- regobj
      outvar.class <- NULL
  }
  class(regs) <- c(class(regs), "ListOfRegressionForms")
  if (censoring) regs <- stratify_by_uncensored(regs)
  return(regs)
}

## --------------------------------------------------------
## DEFINES THE LOWEST LEVEL REPRESENTATION FOR A SINGLE REGRESSION FORMULA BASED WITH THE FOLLOWING STRUCTURE
## --------------------------------------------------------
## sep_predvars_sets = FALSE
## str(all_outVar_nms[[1]])
## lapply(all_outVar_nms[[1]], function(var) {return(var)})
## chr [1:3] "C" "TI.t" "N.t"
## str(all_predVar_nms[[1]])
## chr [1:2] "highA1c.t" "lastN.t"
## str(all_outVar_class[[1]])
## List of 3
## $ C   : chr "binary"
## $ TI.t: chr "binary"
## $ N.t : chr "binary"
## str(all_subsets_vars[[1]])
## List of 3
## $ : chr "C"
## $ : chr "TI.t"
## $ : chr "N.t"
## str(all_subset_exprs[[1]])
## $ : chr "rep.int(TRUE, .N)"
## $ : chr "rep.int(TRUE, .N)"
## $ : chr "rep.int(TRUE, .N)"
SingleRegressionFormClass <- R6Class("SingleRegressionFormClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    outvar = character(),          # vector of regression outcome variable names
    predvars = character(),        # vector of predictor names
    outvar.class = list(),         # Named LIST of outcome class names: binary
    subset_vars = list(),          # Named LIST for subset vars, one list item per outcome in outvar, each list item can be a character vector.
                                   # Later these are tested for missing values, which forms the basis of the logical subset vector)
    subset_exprs = list(),         # Named LIST of subset expressions (as strings), one list item per outcome in outvar.
                                   # Each item is a vector of different subsetting expressions (form stratified models)
                                   # These expressions are evaluated in the envir of the data, must evaluate to a logical vector
    model_contrl = list(),
    censoring = FALSE,             #
    initialize = function(outvar,
                          predvars,
                          outvar.class,
                          subset_vars = NULL,
                          subset_exprs = NULL,
                          model_contrl = NULL,
                          censoring = FALSE) {

      assert_that(is.character(outvar) || is.null(outvar))
      assert_that(is.character(predvars) || is.null(predvars))

      if (is.null(outvar)) {
        outvar <- "NULL"
        outvar.class <- list("NULL")
      }

      self$outvar <- outvar
      self$predvars <- predvars
      self$outvar.class <- self$checkInputList(outvar.class)

      if (!is.null(subset_vars)) {
        self$subset_vars <- self$checkInputList(subset_vars)
      } else {
        self$subset_vars <- lapply(self$outvar, function(var) {var})
        names(self$subset_vars) <- self$outvar
      }

      if (!is.null(subset_exprs)) {
        self$subset_exprs <- self$checkInputList(subset_exprs)
      } else {
        self$subset_exprs <- lapply(self$outvar, function(var) {NULL})
        names(self$subset_exprs) <- self$outvar
      }

      assert_that(is.list(model_contrl))
      if (length(model_contrl)>0) {
        if (any(is.null(names(model_contrl))) || any(names(model_contrl) %in% "")) stop("all items in list 'model_contrl' must be named")
      }

      self$model_contrl <- model_contrl

      assert_that(is.flag(censoring))
      self$censoring <- censoring
      return(self)
    },

    show = function() {
      str(self$get.reg)
      return(invisible(self$get.reg))
    },

    checkInputList = function(inputlist) {
      assert_that(is.list(inputlist))
      assert_that(length(inputlist) == length(self$outvar))
      assert_that(all(names(inputlist) %in% self$outvar))
      return(inputlist)
    }
  ),
  active = list(
    get.reg = function() {
      list(outvar = self$outvar,
           predvars = self$predvars,
           outvar.class = self$outvar.class,
           subset_vars = self$subset_vars,
           subset_exprs = self$subset_exprs,
           model_contrl = self$model_contrl,
           censoring = self$censoring
           )
    }
  )
)

## --------------------------------------------------------
## GENERAL RegressionClass THAT INHERITS FROM SingleRegressionFormClass
## THIS CLASS IS used (and subsetted, if needed) BY all classes that inherit from ModelGeneric class
## --------------------------------------------------------
RegressionClass <- R6Class("RegressionClass",
  inherit = SingleRegressionFormClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard)
    ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    initialize = function(ReplMisVal0 = TRUE, ...) {
      self$ReplMisVal0 <- ReplMisVal0
      super$initialize(...)
    },

    # Take the clone of a parent RegressionClass (reg) for length(self$outvar) regressions
    # and set self to a single univariate k_i regression for outcome self$outvar[[k_i]]
    ChangeManyToOneRegresssion = function(k_i, reg) {
      assert_that(!missing(k_i))
      if (missing(reg)) stop("reg must be also specified when k_i is specified")
      assert_that(is.count(k_i))
      n_regs <- get_n_regs(reg)
      assert_that(k_i <= n_regs)
      # S3 method dispatch on class of reg$RegressionForms:
      select_reg(reg, k_i = k_i, self)
      # self$resetS3class()
      # Setting the self class for S3 dispatch on SummaryModel type
      if ("ListOfRegressionForms" %in% class(reg$RegressionForms)) {
        self$S3class <- "generic" # Multivariate/Univariate regression at the top level, need to do another round of S3 dispatch on SummaryModel
      } else if (length(self$outvar)==1L && length(self$subset_exprs) > 1L) {
        self$S3class <- "stratify" # Set class on outvar.class for S3 dispatch...
      } else if (length(self$outvar)==1L && length(self$subset_exprs) <= 1L) {
        self$S3class <- self$outvar.class # Set class on outvar.class for S3 dispatch...
      } else if (length(self$outvar) > 1){
        stop("can't define a univariate regression for an outcome of length > 1")
      } else {
        stop("can't have an outcome with no class type")
      }
      return(invisible(self))
    }
    # resetS3class = function() class(self) <- c("RegressionClass", "R6")
  ),

  active = list(
    # For S3 dispatch on newsummarymodel():
    S3class = function(newclass) {
      if (!missing(newclass)) {
        # if (length(class(self)) > 2) stop("S3 dispatch class on RegressionClass has already been set")
        if (length(newclass) > 1) stop("cannot set S3 method dispatch to more than 1 class from: generic, contin, categor, stratify, binary or Qlearn")
        private$.S3class <- newclass
        return(invisible(NULL))
      } else {
        S3class <- private$.S3class
        class(S3class) <- S3class
        return(S3class)
      }
    }
  ),
  private = list(
    .S3class = "generic"
  )
)

## --------------------------------------------------------
## RegressionClass that controls a single model (one iteration) in modeling of the iterative G-COMP.
## INHERITS FROM RegressionClass.
## --------------------------------------------------------
RegressionClassQlearn <- R6Class("RegressionClassQlearn",
  inherit = RegressionClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    Qreg_counter = integer(),
    all_Qregs_indx = integer(),
    t_period = integer(),
    TMLE = FALSE,
    CVTMLE = FALSE,
    byfold_Q = FALSE,
    keep_idx = FALSE,           ## should ModelQlearn remove internally stored subset of indices used for training?
    keep_model_fit = TRUE,      ## keep the model fit object for current Q_k
    stratifyQ_by_rule = FALSE,  ## train only among those who are following the rule of interest?
    lower_bound_zero_Q = TRUE,
    skip_update_zero_Q = TRUE,
    regimen_names = NA,
    pool_regimes = FALSE,
    maxpY = 1.0,                ## max incidence P(Y=1|...) for rare-outcomes TMLE, only works with learners that can handle logistic-link with outcomes > 1
    TMLE_updater = "TMLE.updater.speedglm",
    initialize = function(Qreg_counter,
                          all_Qregs_indx,
                          t_period,
                          TMLE,
                          CVTMLE,
                          byfold_Q,
                          stratifyQ_by_rule,
                          regimen_names,
                          pool_regimes,
                          keep_idx,
                          keep_model_fit,
                          lower_bound_zero_Q = stremrOptions("lower_bound_zero_Q"),
                          skip_update_zero_Q = stremrOptions("skip_update_zero_Q"),
                          maxpY,
                          TMLE_updater,
                          ...) {
      self$Qreg_counter <- Qreg_counter
      self$all_Qregs_indx <- all_Qregs_indx
      self$t_period <- t_period

      if (!missing(TMLE)) self$TMLE <- TMLE
      if (!missing(CVTMLE)) self$CVTMLE <- CVTMLE
      if (!missing(byfold_Q)) self$byfold_Q <- byfold_Q
      if (!missing(stratifyQ_by_rule)) self$stratifyQ_by_rule <- stratifyQ_by_rule
      if (!missing(regimen_names)) self$regimen_names <- regimen_names
      if (!missing(pool_regimes)) self$pool_regimes <- pool_regimes
      if (!missing(keep_idx)) self$keep_idx <- keep_idx
      if (!missing(keep_model_fit)) self$keep_model_fit <- keep_model_fit
      if (!missing(maxpY)) self$maxpY <- maxpY
      if (!missing(TMLE_updater)) self$TMLE_updater <- TMLE_updater

      self$lower_bound_zero_Q <- lower_bound_zero_Q
      self$skip_update_zero_Q <- skip_update_zero_Q

      super$initialize(...)
    }
  ),
  active = list(
    get.reg = function() {
      list(Qreg_counter = self$Qreg_counter,
           t_period = self$t_period,
           TMLE = self$TMLE,
           CVTMLE = self$CVTMLE,
           byfold_Q = self$byfold_Q,
           outvar = self$outvar,
           predvars = self$predvars,
           outvar.class = self$outvar.class,
           subset_vars = self$subset_vars,
           subset_exprs = self$subset_exprs,
           subset_censored = self$subset_censored,
           stratifyQ_by_rule = self$stratifyQ_by_rule,
           lower_bound_zero_Q = self$lower_bound_zero_Q,
           regimen_names = self$regimen_names,
           pool_regimes = self$pool_regimes,
           keep_idx = self$keep_idx,
           keep_model_fit = self$keep_model_fit,
           model_contrl = self$model_contrl,
           censoring = self$censoring,
           maxpY = self$maxpY,
           TMLE_updater = self$TMLE_updater
           )
    }
  )
)

RegressionClassSDR <- R6Class("RegressionClassSDR",
  inherit = RegressionClassQlearn,
  public = list(
    SDR_model = NULL,
    stabilize = FALSE,
    initialize = function(SDR_model, stabilize, ...) {
      self$SDR_model <- SDR_model
      self$stabilize <- stabilize
      super$initialize(...)
    }
  )
)

# ---------------------------------------------------------------------------------
# S3 methods for regression subsetting/stratification/interval subsetting
# ---------------------------------------------------------------------------------
select_reg <- function(RegressionForms, reg, k_i, self) { UseMethod("select_reg") }
select_reg.ListOfRegressionForms <- function(RegressionForms, reg, k_i, self) {
  self$RegressionForms <- RegressionForms[[k_i]]
  return(invisible(self))
}

select_reg.RegressionClass <- function(RegressionClassObj, k_i, self) {
  n_regs <- get_n_regs(RegressionClassObj)
  self$outvar.class <- RegressionClassObj$outvar.class[[k_i]]
  self$outvar <- RegressionClassObj$outvar[[k_i]] # An outcome variable that is being modeled
  self$model_contrl <- RegressionClassObj$model_contrl

  if (self$reg_hazard) {
    # Modeling hazards of bin indicators, no need to condition on previous outcomes as they will all be degenerate. P(A,B,C|D) -> P(A|D),P(B|D),P(C|D)
    self$predvars <- RegressionClassObj$predvars # Predictors
  } else {
    # Factorizating the joint prob as P(A,B,C|D):=P(A|D)*P(B|A,D)*P(C|A,B,D)
    self$predvars <- c(RegressionClassObj$outvar[-c(k_i:n_regs)], RegressionClassObj$predvars) # Predictors
    # The subset_vars is a list when RegressionClass is used to specify several regression models.
    # Obtain appropriate subset_vars for this regression (k_i) and set it to self.
    # On the other hand, if subset_vars is a vector of variable names, all of those variables will be used for
    # choosing the subset_vars for all n_regs regressions.
  }

  if (is.list(RegressionClassObj$subset_vars)) {
    self$subset_vars <- RegressionClassObj$subset_vars[[k_i]]
  }

  # Doing the same for subset_exprs for stratified models on subsets
  if (is.list(RegressionClassObj$subset_exprs)) {
    self$subset_exprs <- RegressionClassObj$subset_exprs[[k_i]]
  }

  # Doing the same for intervals when modeling continuous outcomes
  # if (("contin" %in% class(self)) && is.list(reg$intrvls)) {
  # if (is.list(reg$intrvls)) {
  #   outvar_idx <- which(names(reg$intrvls) %in% self$outvar)
  #   self$intrvls <- reg$intrvls[[outvar_idx]]
  # }

  return(invisible(self))
}

# ---------------------------------------------------------------------------
# S3 methods for SingleRegressionFormClass and ListOfRegressionForms
# ---------------------------------------------------------------------------
print.SingleRegressionFormClass <- function(singleregobj) singleregobj$show()

get_n_regs <- function(singleregobj) { UseMethod("get_n_regs") }
get_n_regs.ListOfRegressionForms <- function(regobjlist) return(length(regobjlist))
get_n_regs.SingleRegressionFormClass <- function(singleregobj) return(length(singleregobj$outvar))
# get_n_regs.RegressionClass <- function(regobj) return(get_n_regs(regobj$RegressionForms))
get_n_regs.RegressionClass <- function(regobj) return(length(regobj$outvar))

get_outvars <- function(regobjlist) { UseMethod("get_outvars") }
get_outvars.ListOfRegressionForms <- function(regobjlist) {
  outvars <- NULL
  for (idx in seq_along(regobjlist))
    outvars <- c(outvars, regobjlist[[idx]]$outvar)
  return(outvars)
}
get_subset_exprs <- function(regobjlist) { UseMethod("get_subset_exprs") }
get_subset_exprs.ListOfRegressionForms <- function(regobjlist) {
  subset_exprs <- NULL
  for (idx in seq_along(regobjlist))
    subset_exprs <- c(subset_exprs, regobjlist[[idx]]$subset_exprs)
  return(subset_exprs)
}
set_subset_exprs <- function(regobjlist, idx, subset_exprs) { UseMethod("set_subset_exprs") }
set_subset_exprs.ListOfRegressionForms <- function(regobjlist, idx, subset_exprs) {
  # subset_exprs <- NULL
  idx_count <- 0
  for (idx_reg in seq_along(regobjlist)) {
    for (idx_outvar in seq_along(regobjlist[[idx_reg]]$outvar)) {
      idx_count <- idx_count + 1
      if (idx == idx_count) {
        regobjlist[[idx_reg]]$subset_exprs[[idx_outvar]] <- subset_exprs
        return(invisible(regobjlist[[idx_reg]]))
      }
    }
  }
}
