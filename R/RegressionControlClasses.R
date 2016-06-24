# --------------------------------------------------------
# DEFINES THE LOWEST LEVEL REPRESENTATION FOR A SINGLE REGRESSION FORMULA BASED WITH THE FOLLOWING STRUCTURE
# --------------------------------------------------------
# sep_predvars_sets = FALSE
# str(all_outVar_nms[[1]])
# lapply(all_outVar_nms[[1]], function(var) {return(var)})
# chr [1:3] "C" "TI.t" "N.t"
# str(all_predVar_nms[[1]])
# chr [1:2] "highA1c.t" "lastN.t"
# str(all_outVar_class[[1]])
# List of 3
# $ C   : chr "binary"
# $ TI.t: chr "binary"
# $ N.t : chr "binary"
# str(all_subsets_vars[[1]])
# List of 3
# $ : chr "C"
# $ : chr "TI.t"
# $ : chr "N.t"
# str(all_subset_exprs[[1]])
# $ : chr "rep.int(TRUE, .N)"
# $ : chr "rep.int(TRUE, .N)"
# $ : chr "rep.int(TRUE, .N)"
# SingleRegressionFormClass <- R6Class("SingleRegressionFormClass",
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     outvar = character(),          # vector of regression outcome variable names
#     predvars = character(),        # vector of predictor names
#     outvar.class = list(),         # Named LIST of outcome class names: binary / continuous / categorical
#     subset_vars = list(),          # Named LIST for subset vars, one list item per outcome in outvar, each list item can be a character vector.
#                                    # Later these are tested for missing values, which forms the basis of the logical subset vector)
#     subset_exprs = list(),         # Named LIST of subset expressions (as strings), one list item per outcome in outvar.
#                                    # Each item is a vector of different subsetting expressions (form stratified models)
#                                    # These expressions are evaluated in the envir of the data, must evaluate to a logical vector
#     model_contrl = list(),
#     censoring = FALSE,             #
#     initialize = function(outvar, predvars, outvar.class, subset_vars = NULL, subset_exprs = NULL, model_contrl = NULL, censoring = FALSE) {
#       assert_that(is.character(outvar))
#       assert_that(is.character(predvars) || is.null(predvars))
#       self$outvar <- outvar
#       self$predvars <- predvars
#       self$outvar.class <- self$checkInputList(outvar.class)

#       if (!is.null(subset_vars)) {
#         self$subset_vars <- self$checkInputList(subset_vars)
#       } else {
#         self$subset_vars <- lapply(self$outvar, function(var) {var})
#         names(self$subset_vars) <- self$outvar
#       }

#       if (!is.null(subset_exprs)) {
#         self$subset_exprs <- self$checkInputList(subset_exprs)
#       } else {
#         self$subset_exprs <- lapply(self$outvar, function(var) {NULL})
#         names(self$subset_exprs) <- self$outvar
#       }

#       assert_that(is.list(model_contrl))
#       self$model_contrl <- model_contrl

#       assert_that(is.flag(censoring))
#       self$censoring <- censoring
#       return(self)
#     },

#     show = function() {
#       str(self$get.reg)
#       return(invisible(self$get.reg))
#     },

#     checkInputList = function(inputlist) {
#       assert_that(is.list(inputlist))
#       assert_that(length(inputlist) == length(self$outvar))
#       assert_that(all(names(inputlist) %in% self$outvar))
#       return(inputlist)
#     }
#   ),
#   active = list(
#     get.reg = function() {
#       list(outvar = self$outvar,
#            predvars = self$predvars,
#            outvar.class = self$outvar.class,
#            subset_vars = self$subset_vars,
#            subset_exprs = self$subset_exprs,
#            model_contrl = self$model_contrl,
#            censoring = self$censoring
#            )
#     }
#   )
# )

# RegressionClass <- R6Class("RegressionClass",
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     RegressionForms = NULL,
#     reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard)
#     ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
#     nbins = NULL,                  # actual nbins used, for cont. outvar, defined in ContinModel$new()
#     bin_nms = NULL,                # column names for bin indicators

#     # useglm = logical(),            # TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
#     # GLMpackage = character(),      # 'glm', 'speedglm' or 'h2o'
#     # fit.package = character(),
#     # fit.algorithm = character(),
#     fit.package = c("speedglm", "glm", "h2o"),
#     fit.algorithm = c("GLM", "GBM", "RF", "SL"),

#     parfit = logical(),            # TRUE for fitting binary regressions in parallel
#     bin_bymass = logical(),        # for cont outvar, create bin cutoffs based on equal mass distribution?
#     bin_bydhist = logical(),       # if TRUE, use dhist approach for bin definitions
#     max_nperbin = integer(),       # maximum n observations allowed per binary bin
#     pool_cont = logical(),         # Pool binned cont outvar obs into long format (adding bin_ID as a covaraite)
#     outvars_to_pool = character(), # Names of the binned continuous sVars, should match bin_nms
#     intrvls.width = 1L,            # Named vector of bin-widths (bw_j : j=1,...,M) for each each bin in self$intrvls
#                                    # When sA is not continuous, intrvls.width IS SET TO 1.
#                                    # When sA is continuous, intrvls.width is SET TO self$intrvls.width INSIDE ContinModel$new() with names(intrvls.width) <- reg$bin_nms
#                                    # CAN BE QUERIED BY BinaryOutcomeModel$predictAeqa() as: intrvls.width[outvar]
#     intrvls = numeric(),           # Vector of numeric cutoffs defining the bins or a named list of numeric intervals (for length(self$outvar) > 1)
#     levels = NULL,
#     # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
#     initialize = function(RegressionForms,
#                           intrvls,
#                           ReplMisVal0 = TRUE, # Needed to add ReplMisVal0 = TRUE for case sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise:
#                           fit.package = getopt("fit.package"),
#                           fit.algorithm = getopt("fit.algorithm"),
#                           parfit = getopt("parfit"),
#                           nbins = getopt("nbins"),
#                           bin_bymass = getopt("bin.method")%in%"equal.mass",
#                           bin_bydhist = getopt("bin.method")%in%"dhist",
#                           max_nperbin = getopt("maxNperBin"),
#                           pool_cont = getopt("poolContinVar")
#                           ) {

#       if (!missing(RegressionForms)) self$RegressionForms <- RegressionForms
#       self$ReplMisVal0 <- ReplMisVal0
#       self$fit.package <- fit.package
#       self$fit.algorithm <- fit.algorithm

#       self$parfit <- parfit
#       self$nbins <- nbins
#       self$bin_bymass <- bin_bymass
#       self$bin_bydhist <- bin_bydhist
#       self$max_nperbin <- max_nperbin
#       self$pool_cont <- pool_cont

#       if (!missing(intrvls)) {
#         assert_that(is.list(intrvls))
#         assert_that(length(outvar) == length(intrvls))
#         assert_that(all(names(intrvls) %in% outvar))
#         self$intrvls <- intrvls
#       } else {
#         self$intrvls <- NULL
#       }
#       self$levels <- NULL
#     },

#     # take the clone of a parent RegressionClass (reg) for length(self$outvar) regressions
#     # and set self to a single univariate k_i regression for outcome self$outvar[[k_i]]
#     ChangeManyToOneRegresssion = function(k_i, reg) {
#       assert_that(!missing(k_i))
#       if (missing(reg)) stop("reg must be also specified when k_i is specified")
#       assert_that(is.count(k_i))
#       n_regs <- get_n_regs(reg)
#       assert_that(k_i <= n_regs)
#       # S3 method dispatch on class of reg$RegressionForms:
#       select_reg(reg$RegressionForms, reg = reg, k_i = k_i, self)
#       self$resetS3class()
#       # Setting the self class for S3 dispatch on SummaryModel type
#       if ("ListOfRegressionForms" %in% class(reg$RegressionForms)) {
#         self$S3class <- "generic" # Multivariate/Univariate regression at the top level, need to do another round of S3 dispatch on SummaryModel
#       } else if (length(self$outvar)==1L && length(self$subset_exprs) > 1L) {
#         self$S3class <- "stratify" # Set class on outvar.class for S3 dispatch...
#       } else if (length(self$outvar)==1L && length(self$subset_exprs) <= 1L) {
#         self$S3class <- self$outvar.class # Set class on outvar.class for S3 dispatch...
#       } else if (length(self$outvar) > 1){
#         stop("can't define a univariate regression for an outcome of length > 1")
#       } else {
#         stop("can't have an outcome with no class type")
#       }
#       return(invisible(self))
#     },
#     show = function(print_format = TRUE) {
#       print(str(self$get.reg))
#       return(invisible(self))
#     },
#     resetS3class = function() class(self) <- c("RegressionClass", "R6")
#   ),

#   active = list(
#     # For S3 dispatch on newsummarymodel():
#     S3class = function(newclass) {
#       if (!missing(newclass)) {
#         if (length(class(self)) > 2) stop("S3 dispatch class on RegressionClass has already been set")
#         if (length(newclass) > 1) stop("cannot set S3 class on RegressionClass with more than one outvar variable")
#         class(self) <- c(class(self), newclass)
#       } else {
#         return(class(self))
#       }
#     },
#     outvar.class = function(outvar.class) {
#       if (missing(outvar.class)) {
#         return(self$RegressionForms$outvar.class)
#       } else {
#         self$RegressionForms$outvar.class <- outvar.class
#         return(invisible(self))
#       }
#     },
#     outvar = function(outvar) {
#       if (missing(outvar)) {
#         return(self$RegressionForms$outvar)
#       } else {
#         self$RegressionForms$outvar <- outvar
#         return(invisible(self))
#       }
#     },
#     predvars = function(predvars) {
#       if (missing(predvars)) {
#         return(self$RegressionForms$predvars)
#       } else {
#         self$RegressionForms$predvars <- predvars
#         return(invisible(self))
#       }
#     },
#     subset_vars = function(subset_vars) {
#       if (missing(subset_vars)) {
#         return(self$RegressionForms$subset_vars)
#       } else {
#         self$RegressionForms$subset_vars <- subset_vars
#         return(invisible(self))
#       }
#     },
#     subset_exprs = function(subset_exprs) {
#       if (missing(subset_exprs)) {
#         return(self$RegressionForms$subset_exprs)
#       } else {
#         self$RegressionForms$subset_exprs <- subset_exprs
#         return(invisible(self))
#       }
#     },
#     model_contrl = function(model_contrl) {
#       if (missing(model_contrl)) {
#         return(self$RegressionForms$model_contrl)
#       } else {
#         self$RegressionForms$model_contrl <- model_contrl
#         return(invisible(self))
#       }
#     },
#     get.reg = function() { self$RegressionForms$get.reg }
#   )
# )


# ---------------------------------------------------------------------------------
# S3 methods for regression subsetting/stratification/interval subsetting
# ---------------------------------------------------------------------------------
select_reg <- function(RegressionForms, reg, k_i, self) { UseMethod("select_reg") }
select_reg.ListOfRegressionForms <- function(RegressionForms, reg, k_i, self) {
  self$RegressionForms <- RegressionForms[[k_i]]
  return(invisible(self))
}

# select_reg.SingleRegressionFormClass <- function(RegressionClassObj, reg, k_i, self) {
#   n_regs <- get_n_regs(RegressionClassObj)
#   self$RegressionClassObj <- RegressionClassObj$clone()
#   self$outvar.class <- RegressionClassObj$outvar.class[[k_i]]
#   self$outvar <- RegressionClassObj$outvar[[k_i]] # An outcome variable that is being modeled
#   self$model_contrl <- RegressionClassObj$model_contrl

#   if (self$reg_hazard) {
#     # Modeling hazards of bin indicators, no need to condition on previous outcomes as they will all be degenerate. P(A,B,C|D) -> P(A|D),P(B|D),P(C|D)
#     self$predvars <- RegressionClassObj$predvars # Predictors
#   } else {
#     # Factorizating the joint prob as P(A,B,C|D):=P(A|D)*P(B|A,D)*P(C|A,B,D)
#     self$predvars <- c(RegressionClassObj$outvar[-c(k_i:n_regs)], RegressionClassObj$predvars) # Predictors
#     # The subset_vars is a list when RegressionClass is used to specify several regression models.
#     # Obtain appropriate subset_vars for this regression (k_i) and set it to self.
#     # On the other hand, if subset_vars is a vector of variable names, all of those variables will be used for
#     # choosing the subset_vars for all n_regs regressions.
#   }
#   if (is.list(RegressionClassObj$subset_vars)) {
#     self$subset_vars <- RegressionClassObj$subset_vars[[k_i]]
#   }
#   # Doing the same for subset_exprs for stratified models on subsets
#   if (is.list(RegressionClassObj$subset_exprs)) {
#     self$subset_exprs <- RegressionClassObj$subset_exprs[[k_i]]
#   }
#   # Doing the same for intervals when modeling continuous outcomes
#   # if (("contin" %in% class(self)) && is.list(reg$intrvls)) {
#   if (is.list(reg$intrvls)) {
#     outvar_idx <- which(names(reg$intrvls) %in% self$outvar)
#     self$intrvls <- reg$intrvls[[outvar_idx]]
#   }
#   return(invisible(self))
# }

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
  subset_exprs <- NULL
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



