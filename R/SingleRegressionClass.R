
# ---------------------------------------------------------------------------
# S3 methods for SingleRegressionFormClass and ListOfRegressionForms
# ---------------------------------------------------------------------------
print.SingleRegressionFormClass <- function(singleregobj) singleregobj$show()

get_n_regs <- function(singleregobj) { UseMethod("get_n_regs") }
get_n_regs.ListOfRegressionForms <- function(regobjlist) return(length(regobjlist))
get_n_regs.SingleRegressionFormClass <- function(singleregobj) return(length(singleregobj$outvar))
get_n_regs.RegressionClass <- function(regobj) return(get_n_regs(regobj$RegressionForms))

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
set_subset_exprs <- function(regobjlist, idx, subset_expr) { UseMethod("set_subset_exprs") }
set_subset_exprs.ListOfRegressionForms <- function(regobjlist, idx, subset_expr) {
  subset_exprs <- NULL
  idx_count <- 0
  for (idx_reg in seq_along(regobjlist)) {
    for (idx_outvar in seq_along(regobjlist[[idx_reg]]$outvar)) {
      idx_count <- idx_count + 1
      if (idx == idx_count) {
        regobjlist[[idx_reg]]$subset_exprs[[idx_outvar]] <- subset_expr
        return(invisible(regobjlist[[idx_reg]]))
      }
    }
  }
}

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
SingleRegressionFormClass <- R6Class("SingleRegressionFormClass",
  class = TRUE,
  # inherit = StratifiedRegressionModelClass,
  portable = TRUE,
  public = list(
    outvar = character(),          # vector of regression outcome variable names
    predvars = character(),        # vector of predictor names
    outvar.class = list(),         # Named LIST of outcome class names: binary / continuous / categorical
    subset_vars = list(),          # Named LIST for subset vars, one list item per outcome in outvar, each list item can be a character vector.
                                   # Later these are tested for missing values, which forms the basis of the logical subset vector)
    subset_exprs = list(),         # Named LIST of subset expressions (as strings), one list item per outcome in outvar.
                                   # Each item is a vector of different subsetting expressions (form stratified models)
                                   # These expressions are evaluated in the envir of the data, must evaluate to a logical vector
    censoring = FALSE,             #
    formula = NULL,                # (NOT IMPLEMENTED) reg formula, if provided run using the usual glm / speedglm functions
    # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
    initialize = function(outvar, predvars, outvar.class, subset_vars = NULL, subset_exprs = NULL, censoring = FALSE, formula = NULL) {
      assert_that(is.character(outvar))
      assert_that(is.character(predvars) || is.null(predvars))
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

      assert_that(is.flag(censoring))
      self$censoring <- censoring
      if (!is.null(formula)) self$formula <- formula
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
           censoring = self$censoring
           )
    }
    # get.reg = function() {
    #   c(super$get.reg, censoring = self$censoring)
  )
)



# ContinRegressionFormClass <- R6Class("ContinRegressionFormClass",
#   class = TRUE,
#   inherit = SingleRegressionFormClass,
#   portable = TRUE
# )

# CatRegressionFormClass <- R6Class("ContinRegressionFormClass",
#   class = TRUE,
#   inherit = SingleRegressionFormClass,
#   portable = TRUE
# )

# BinRegressionFormClass <- R6Class("BinRegressionFormClass",
#   class = TRUE,
#   inherit = SingleRegressionFormClass,
#   portable = TRUE
# )
# StratifiedRegressionModelClass <- R6Class("StratifiedRegressionModelClass",
#   class = TRUE,
#   portable = TRUE,
#   public = list(
#     outvar = character(),          # Single outcome variable name
#     predvars = character(),        # Vector of predictor names
#     outvar.class = character(),    # Character class name for outcome: binary / continuous / categorical
#     subset_vars = character(),     # Character vector for subsettting var names. Later these are tested for missing values, which forms the basis of the logical subset vector)
#     subset_exprs = character(),    # Character string for a single stratification rule (this rule defines a new model for subsetting)
#                                    # This expression will be evaluated in the envir of the data and must evaluate to a logical vector
#     formula = NULL,                # (NOT IMPLEMENTED) reg formula, if provided run using the usual glm / speedglm functions
#     # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
#     initialize = function(outvar, predvars, outvar.class, subset_vars = NULL, subset_exprs = NULL, formula = NULL) {
#       assert_that(is.character(outvar))
#       assert_that(length(outvar)==1)
#       assert_that(is.character(predvars) || is.null(predvars))
#       assert_that(is.character(outvar.class))
#       assert_that(length(outvar.class==1L))

#       self$outvar <- outvar
#       self$predvars <- predvars
#       self$outvar.class <- outvar.class

#       if (!is.null(subset_vars)) {
#         assert_that(is.character(subset_vars))
#         self$subset_vars <- subset_vars
#       } else {
#         self$subset_vars <- NULL
#       }
#       if (!is.null(subset_exprs)) {
#         assert_that(is.character(subset_exprs))
#         assert_that(length(subset_exprs)==1L)
#         self$subset_exprs <- subset_exprs
#       } else {
#         self$subset_exprs <- NULL
#       }
#       if (!is.null(formula)) self$formula <- formula
#       return(self)
#     },

#     show = function() {
#       str(self$get.reg)
#       return(invisible(self$get.reg))
#     }
#   ),
#   active = list(
#     get.reg = function() {
#       list(outvar = self$outvar,
#            predvars = self$predvars,
#            outvar.class = self$outvar.class,
#            subset_vars = self$subset_vars,
#            subset_exprs = self$subset_exprs
#            )
#     }
#   )
# )