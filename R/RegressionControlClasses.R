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
    model_contrl = list(),
    censoring = FALSE,             #
    initialize = function(outvar, predvars, outvar.class, subset_vars = NULL, subset_exprs = NULL, model_contrl = NULL, censoring = FALSE) {
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

# --------------------------------------------------------
# MORE GENERAL RegressionClass THAT INHERITS FROM SingleRegressionFormClass AND IS CLOSED AND SUBSETTED BY GENERIC MODEL
# --------------------------------------------------------

## ---------------------------------------------------------------------
#' R6 class that defines regression models evaluating P(sA|sW), for summary measures (sW,sA)
#'
#' This R6 class defines fields and methods that controls all the parameters for non-parametric
#'  modeling and estimation of multivariate joint conditional probability model \code{P(sA|sW)} for summary measures \code{(sA,sW)}.
#'  Note that \code{sA} can be multivariate and any component of \code{sA[j]} can be either binary, categorical or continuous.
#'  The joint probability for \code{P(sA|sA)} = \code{P(sA[1],...,sA[k]|sA)} is first factorized as
#'  \code{P(sA[1]|sA)} * \code{P(sA[2]|sA, sA[1])} * ... * \code{P(sA[k]|sA, sA[1],...,sA[k-1])},
#'  where each of these conditional probability models is defined by a new instance of a \code{\link{GenericModel}} class
#'  (and a corresponding instance of the \code{RegressionClass} class).
#'  If \code{sA[j]} is binary, the conditional probability \code{P(sA[j]|sW,sA[1],...,sA[j-1])} is evaluated via logistic regression model.
#'  When \code{sA[j]} is continuous (or categorical), its estimation will be controlled by a new instance of
#'  the \code{\link{ContinModel}} class (or the \code{\link{CategorModel}} class), as well as the accompanying new instance of the
#'  \code{RegressionClass} class. The range of continuous \code{sA[j]} will be fist partitioned into \code{K} bins and the corresponding \code{K}
#'  bin indicators (\code{B_1,...,B_K}), with \code{K} new instances of \code{\link{GenericModel}} class, each instance defining a
#'  single logistic regression model for one binary bin indicator outcome \code{B_j} and predictors (\code{sW, sA[1],...,sA[k-1]}).
#'  Thus, the first instance of \code{RegressionClass} and \code{GenericModel} classes will automatically
#'  spawn recursive calls to new instances of these classes until the entire tree of binary logistic regressions that defines
#'  the joint probability \code{P(sA|sW)} is build.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{sep_predvars_sets}} - Logical indicating the type of regression to run,
#'    if \code{TRUE} fit the joint P(\code{outvar}|\code{predvars}) (default),
# '   if \code{FALSE}, fit P(\code{outvar[1]}|\code{predvars[[1]]})*...*P(\code{outvar[K]}|\code{predvars[[K]}].
#'    More specifically, if \code{FALSE} (default), use the same predictors in \code{predvars} (vector of names) for all nodes in \code{outvar};
#'    when \code{TRUE} uses separate sets in \code{predvars} (must be a named list of character vectors) for fitting each node in \code{outvar}.
#' \item{\code{outvar.class}} - Character vector indicating a class of each outcome var: \code{bin} / \code{cont} / \code{cat}.
#' \item{\code{outvar}} - Character vector of regression outcome variable names.
#' \item{\code{predvars}} - Either a pooled character vector of all predictors (\code{sW}) or a vector of regression-specific predictor names.
#'      When \code{sep_predvars_sets=TRUE}, this must be a named list of predictor names, the list names corresponding to each node name in \code{outvar},
#'      and each list item being a vector specifying the regression predictors for a specific outcome in \code{outvar}.
#' \item{{reg_hazard}} - Logical, if TRUE, the joint probability model P(outvar | predvars) is factorized as
#'    \\prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard).
#' \item{\code{subset_vars}} - Subset variables (later evaluated to logical vector based on non-missing (!is.na()) values of these variables).
#' \item{\code{subset_exprs}} - Subset expressions (later evaluated to logical vector in the envir of the data).
#' \item{\code{ReplMisVal0}} - Logical, if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0).
#' \item{\code{nbins}} - Integer number of bins used for a continuous outvar, the intervals are defined inside
#'  \code{ContinModel$new()} and then saved in this field.
#' \item{\code{bin_nms}} - Character vector of column names for bin indicators.
#' \item{\code{useglm}} - Logical, if TRUE then fit the logistic regression model using \code{\link{glm.fit}},
#'    if FALSE use \code{\link{speedglm.wfit}}..
#' \item{\code{parfit}} - Logical, if TRUE then use parallel \code{foreach::foreach} loop to fit and predict binary logistic
#'    regressions (requires registering back-end cluster prior to calling the fit/predict functions)..
#' \item{\code{bin_bymass}} - Logical, for continuous outvar, create bin cutoffs based on equal mass distribution.
#' \item{\code{bin_bydhist}} - Logical, if TRUE, use dhist approach for bin definitions.  See Denby and Mallows "Variations on the
#'    Histogram" (2009)) for more..
#' \item{\code{max_nperbin}} - Integer, maximum number of observations allowed per one bin.
#' \item{\code{pool_cont}} - Logical, pool binned continuous outvar observations across bins and only fit only regression model
#'    across all bins (adding bin_ID as an extra covaraite)..
#' \item{\code{outvars_to_pool}} - Character vector of names of the binned continuous outvars, should match \code{bin_nms}.
#' \item{\code{intrvls.width}} - Named numeric vector of bin-widths (\code{bw_j : j=1,...,M}) for each each bin in \code{self$intrvls}.
#'    When \code{sA} is not continuous, \code{intrvls.width} IS SET TO 1. When sA is continuous and this variable \code{intrvls.width}
#'    is not here, the intervals are determined inside \code{ContinModel$new()} and are assigned to this variable as a list,
#'    with \code{names(intrvls.width) <- reg$bin_nms}. Can be queried by \code{BinaryOutcomeModel$predictAeqa()} as: \code{intrvls.width[outvar]}.
#' \item{\code{intrvls}} - Numeric vector of cutoffs defining the bins or a named list of numeric intervals for \code{length(self$outvar) > 1}.
#' \item{\code{cat.levels}} - Numeric vector of all unique values in categorical outcome variable.
#'    Set by \code{\link{CategorModel}} constructor.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(sep_predvars_sets = FALSE,
#'                   outvar.class = gvars$sVartypes$bin,
#'                   outvar, predvars, subset_vars, subset_exprs, intrvls,
#'                   ReplMisVal0 = TRUE,
#'                   useglm = getopt("useglm"),
#'                   parfit = getopt("parfit"),
#'                   nbins = getopt("nbins"),
#'                   bin_bymass = getopt("bin.method")%in%"equal.mass",
#'                   bin_bydhist = getopt("bin.method")%in%"dhist",
#'                   max_nperbin = getopt("maxNperBin"),
#'                   pool_cont = getopt("poolContinVar")}}{Uses the arguments to instantiate an object of R6 class and define the future regression model.}
#'   \item{\code{ChangeManyToOneRegresssion(k_i, reg)}}{ Take a clone of a parent \code{RegressionClass} (\code{reg}) for \code{length(self$outvar)} regressions
#'    and set self to a single univariate \code{k_i} regression for outcome \code{self$outvar[[k_i]]}.}
#'   \item{\code{ChangeOneToManyRegresssions(regs_list)}}{ Take the clone of a parent \code{RegressionClass} for univariate (continuous outvar) regression
#'     and set self to \code{length(regs_list)} bin indicator outcome regressions.}
#'   \item{\code{resetS3class()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{S3class}}{...}
#'   \item{\code{get.reg}}{...}
#' }
#' @export
RegressionClass <- R6Class("RegressionClass",
  inherit = SingleRegressionFormClass,
  class = TRUE,
  portable = TRUE,
  public = list(
    reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard)
    ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    fit.package = c("speedglm", "glm", "h2o"),
    fit.algorithm = c("glm", "gbm", "randomForest", "SL"),
    parfit = logical(),            # TRUE for fitting binary regressions in parallel
    # Needed to add ReplMisVal0 = TRUE for case sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise:
    initialize = function(ReplMisVal0 = TRUE,
                          fit.package = getopt("fit.package"),
                          fit.algorithm = getopt("fit.algorithm"),
                          parfit = getopt("parfit"), ...) {
      self$ReplMisVal0 <- ReplMisVal0
      self$fit.package <- fit.package
      self$fit.algorithm <- fit.algorithm
      self$parfit <- parfit
      super$initialize(...)
    },

    # take the clone of a parent RegressionClass (reg) for length(self$outvar) regressions
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
