#----------------------------------------------------------------------------------
# Classes that control modelling of the multivariate joint probability model P(sA|sW).
#----------------------------------------------------------------------------------
#' @importFrom assertthat assert_that

# ---------------------------------------------------------------------------------
# S3 constructors for the summary model classes:
# ---------------------------------------------------------------------------------
newsummarymodel <- function(reg, DataStorageClass.g0, ...) { UseMethod("newsummarymodel") }
# Summary model constructor for generic regression with multivariate outcome, but one set of predictors
newsummarymodel.generic <- function(reg, DataStorageClass.g0, ...) SummariesModel$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
# Summary model constructor for binary outcome sA[j]:
newsummarymodel.binary <- function(reg, ...) BinOutModel$new(reg = reg, ...)
# Summary model constructor for continuous outcome sA[j]:
newsummarymodel.contin <- function(reg, DataStorageClass.g0, ...) ContinSummaryModel$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
# Summary model constructor for categorical outcome sA[j]:
newsummarymodel.categor <- function(reg, DataStorageClass.g0, ...) CategorSummaryModel$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)
# Summary model constructor for stratification (by reg$subset_exprs):
newsummarymodel.stratify <- function(reg, DataStorageClass.g0, ...) StratifySummariesModel$new(reg = reg, DataStorageClass.g0 = DataStorageClass.g0, ...)

# ---------------------------------------------------------------------------------
# S3 methods for regression subsetting/stratification/interval subsetting
# ---------------------------------------------------------------------------------
select_reg <- function(RegressionForms, reg, k_i, self) { UseMethod("select_reg") }
select_reg.ListOfRegressionForms <- function(RegressionForms, reg, k_i, self) {
  self$RegressionForms <- RegressionForms[[k_i]]
  return(invisible(self))
}
select_reg.SingleRegressionFormClass <- function(RegressionForms, reg, k_i, self) {
  n_regs <- get_n_regs(RegressionForms)
  self$RegressionForms <- RegressionForms$clone()
  self$outvar.class <- RegressionForms$outvar.class[[k_i]]
  self$outvar <- RegressionForms$outvar[[k_i]] # An outcome variable that is being modeled
  if (self$reg_hazard) {
    # Modeling hazards of bin indicators, no need to condition on previous outcomes as they will all be degenerate. P(A,B,C|D) -> P(A|D),P(B|D),P(C|D)
    self$predvars <- RegressionForms$predvars # Predictors
  } else {
    # Factorizating the joint prob as P(A,B,C|D):=P(A|D)*P(B|A,D)*P(C|A,B,D)
    self$predvars <- c(RegressionForms$outvar[-c(k_i:n_regs)], RegressionForms$predvars) # Predictors
    # The subset_vars is a list when RegressionClass is used to specify several regression models.
    # Obtain appropriate subset_vars for this regression (k_i) and set it to self.
    # On the other hand, if subset_vars is a vector of variable names, all of those variables will be used for
    # choosing the subset_vars for all n_regs regressions.
  }
  if (is.list(RegressionForms$subset_vars)) {
    self$subset_vars <- RegressionForms$subset_vars[[k_i]]
  }
  # Doing the same for subset_exprs for stratified models on subsets
  if (is.list(RegressionForms$subset_exprs)) {
    self$subset_exprs <- RegressionForms$subset_exprs[[k_i]]
  }
  # Doing the same for intervals when modeling continuous outcomes
  # if (("contin" %in% class(self)) && is.list(reg$intrvls)) {
  if (is.list(reg$intrvls)) {
    outvar_idx <- which(names(reg$intrvls) %in% self$outvar)
    self$intrvls <- reg$intrvls[[outvar_idx]]
  }
  return(invisible(self))
}

## ---------------------------------------------------------------------
#' R6 class that defines regression models evaluating P(sA|sW), for summary measures (sW,sA)
#'
#' This R6 class defines fields and methods that controls all the parameters for non-parametric
#'  modeling and estimation of multivariate joint conditional probability model \code{P(sA|sW)} for summary measures \code{(sA,sW)}.
#'  Note that \code{sA} can be multivariate and any component of \code{sA[j]} can be either binary, categorical or continuous.
#'  The joint probability for \code{P(sA|sA)} = \code{P(sA[1],...,sA[k]|sA)} is first factorized as
#'  \code{P(sA[1]|sA)} * \code{P(sA[2]|sA, sA[1])} * ... * \code{P(sA[k]|sA, sA[1],...,sA[k-1])},
#'  where each of these conditional probability models is defined by a new instance of a \code{\link{SummariesModel}} class
#'  (and a corresponding instance of the \code{RegressionClass} class).
#'  If \code{sA[j]} is binary, the conditional probability \code{P(sA[j]|sW,sA[1],...,sA[j-1])} is evaluated via logistic regression model.
#'  When \code{sA[j]} is continuous (or categorical), its estimation will be controlled by a new instance of
#'  the \code{\link{ContinSummaryModel}} class (or the \code{\link{CategorSummaryModel}} class), as well as the accompanying new instance of the
#'  \code{RegressionClass} class. The range of continuous \code{sA[j]} will be fist partitioned into \code{K} bins and the corresponding \code{K}
#'  bin indicators (\code{B_1,...,B_K}), with \code{K} new instances of \code{\link{SummariesModel}} class, each instance defining a
#'  single logistic regression model for one binary bin indicator outcome \code{B_j} and predictors (\code{sW, sA[1],...,sA[k-1]}).
#'  Thus, the first instance of \code{RegressionClass} and \code{SummariesModel} classes will automatically
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
#'  \code{ContinSummaryModel$new()} and then saved in this field.
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
#'    is not here, the intervals are determined inside \code{ContinSummaryModel$new()} and are assigned to this variable as a list,
#'    with \code{names(intrvls.width) <- reg$bin_nms}. Can be queried by \code{BinOutModel$predictAeqa()} as: \code{intrvls.width[outvar]}.
#' \item{\code{intrvls}} - Numeric vector of cutoffs defining the bins or a named list of numeric intervals for \code{length(self$outvar) > 1}.
#' \item{\code{cat.levels}} - Numeric vector of all unique values in categorical outcome variable.
#'    Set by \code{\link{CategorSummaryModel}} constructor.
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
  class = TRUE,
  portable = TRUE,
  public = list(
    RegressionForms = NULL,
    reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar (for fitting hazard)
    ReplMisVal0 = TRUE,            # if TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0)
    nbins = NULL,                  # actual nbins used, for cont. outvar, defined in ContinSummaryModel$new()
    bin_nms = NULL,                # column names for bin indicators
    useglm = logical(),            # TRUE to fit reg with glm.fit(), FALSE to fit with speedglm.wfit
    GLMpackage = character(),      # 'glm', 'speedglm' or 'h2o'
    parfit = logical(),            # TRUE for fitting binary regressions in parallel
    bin_bymass = logical(),        # for cont outvar, create bin cutoffs based on equal mass distribution?
    bin_bydhist = logical(),       # if TRUE, use dhist approach for bin definitions
    max_nperbin = integer(),       # maximum n observations allowed per binary bin
    pool_cont = logical(),         # Pool binned cont outvar obs into long format (adding bin_ID as a covaraite)
    outvars_to_pool = character(), # Names of the binned continuous sVars, should match bin_nms
    intrvls.width = 1L,            # Named vector of bin-widths (bw_j : j=1,...,M) for each each bin in self$intrvls
                                   # When sA is not continuous, intrvls.width IS SET TO 1.
                                   # When sA is continuous, intrvls.width is SET TO self$intrvls.width INSIDE ContinSummaryModel$new() with names(intrvls.width) <- reg$bin_nms
                                   # CAN BE QUERIED BY BinOutModel$predictAeqa() as: intrvls.width[outvar]
    intrvls = numeric(),           # Vector of numeric cutoffs defining the bins or a named list of numeric intervals (for length(self$outvar) > 1)
    levels = NULL,
    # family = NULL,               # (NOT IMPLEMENTED) to run w/ other than "binomial" family
    initialize = function(
                          RegressionForms,
                          intrvls,
                          ReplMisVal0 = TRUE, # Needed to add ReplMisVal0 = TRUE for case sA = (netA, sA[j]) with sA[j] continuous, was causing an error otherwise:
                          useglm = getopt("useglm"),
                          GLMpackage = getopt("GLMpackage"),
                          parfit = getopt("parfit"),
                          nbins = getopt("nbins"),
                          bin_bymass = getopt("bin.method")%in%"equal.mass",
                          bin_bydhist = getopt("bin.method")%in%"dhist",
                          max_nperbin = getopt("maxNperBin"),
                          pool_cont = getopt("poolContinVar")
                          ) {

      if (!missing(RegressionForms)) self$RegressionForms <- RegressionForms
      self$ReplMisVal0 <- ReplMisVal0
      self$useglm <- useglm
      self$GLMpackage <- GLMpackage
      self$parfit <- parfit
      self$nbins <- nbins
      self$bin_bymass <- bin_bymass
      self$bin_bydhist <- bin_bydhist
      self$max_nperbin <- max_nperbin
      self$pool_cont <- pool_cont

      if (!missing(intrvls)) {
        assert_that(is.list(intrvls))
        assert_that(length(outvar) == length(intrvls))
        assert_that(all(names(intrvls) %in% outvar))
        self$intrvls <- intrvls
      } else {
        self$intrvls <- NULL
      }
      self$levels <- NULL
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
      select_reg(reg$RegressionForms, reg = reg, k_i = k_i, self)
      self$resetS3class()
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
    },
    show = function(print_format = TRUE) {
      print(str(self$get.reg))
      return(invisible(self))
    },
    resetS3class = function() class(self) <- c("RegressionClass", "R6")
  ),

  active = list(
    # For S3 dispatch on newsummarymodel():
    S3class = function(newclass) {
      if (!missing(newclass)) {
        if (length(class(self)) > 2) stop("S3 dispatch class on RegressionClass has already been set")
        if (length(newclass) > 1) stop("cannot set S3 class on RegressionClass with more than one outvar variable")
        class(self) <- c(class(self), newclass)
      } else {
        return(class(self))
      }
    },
    outvar.class = function(outvar.class) {
      if (missing(outvar.class)) {
        return(self$RegressionForms$outvar.class)
      } else {
        self$RegressionForms$outvar.class <- outvar.class
        return(invisible(self))
      }
    },
    outvar = function(outvar) {
      if (missing(outvar)) {
        return(self$RegressionForms$outvar)
      } else {
        self$RegressionForms$outvar <- outvar
        return(invisible(self))
      }
    },
    predvars = function(predvars) {
      if (missing(predvars)) {
        return(self$RegressionForms$predvars)
      } else {
        self$RegressionForms$predvars <- predvars
        return(invisible(self))
      }
    },
    subset_vars = function(subset_vars) {
      if (missing(subset_vars)) {
        return(self$RegressionForms$subset_vars)
      } else {
        self$RegressionForms$subset_vars <- subset_vars
        return(invisible(self))
      }
    },
    subset_exprs = function(subset_exprs) {
      if (missing(subset_exprs)) {
        return(self$RegressionForms$subset_exprs)
      } else {
        self$RegressionForms$subset_exprs <- subset_exprs
        return(invisible(self))
      }
    },
    formula = function(formula) {
      if (missing(formula)) {
        return(self$RegressionForms$formula)
      } else {
        self$RegressionForms$formula <- formula
        return(invisible(self))
      }
    },
    get.reg = function() { self$RegressionForms$get.reg }
  )
)

## ---------------------------------------------------------------------
#' R6 class for fitting and predicting model P(sA|sW) under g.star or g.0
#'
#' This R6 class Class for defining, fitting and predicting the probability model
#'  \code{P(sA|sW)} under \code{g_star} or under \code{g_0} for summary measures
#'  (\code{sW,sA}). Defines and manages the factorization of the multivariate conditional
#'  probability model \code{P(sA=sa|...)} into univariate regression models
#'  \code{sA[j] ~ sA[j-1] + ... + sA[1] + sW}. The class \code{self$new} method automatically
#'  figures out the correct joint probability factorization into univariate conditional
#'  probabilities based on name ordering provided by (\code{sA_nms}, \code{sW_nms}).
#'  When the outcome variable \code{sA[j]} is binary, this class will automatically call
#'  a new instance of \code{\link{BinOutModel}} class.
#'  Provide \code{self$fit()} function argument \code{data} as a \code{\link{DataStorageClass}} class object.
#'  This data will be used for fitting the model \code{P(sA|sW)}.
#'  Provide \code{self$fit()} function argument \code{newdata} (also as \code{DataStorageClass} class) for predictions of the type
#'  \code{P(sA=1|sW=sw)}, where \code{sw} values are coming from \code{newdata} object.
#'  Finally, provide \code{self$predictAeqa} function \code{newdata} argument
#'  (also \code{DataStorageClass} class object) for getting the likelihood predictions \code{P(sA=sa|sW=sw)}, where
#'  both, \code{sa} and \code{sw} values are coming from \code{newdata} object.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{n_regs}} - .
#' \item{\code{parfit_allowed}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, ...)}}{...}
#'   \item{\code{length}}{...}
#'   \item{\code{getPsAsW.models}}{...}
#'   \item{\code{getcumprodAeqa}}{...}
#'   \item{\code{copy.fit(SummariesModel)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata, ...)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{wipe.alldat}}{...}
#' }
#' @export
SummariesModel <- R6Class(classname = "SummariesModel",
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),      # outcome name(s)
    predvars = character(),    # names of predictor vars
    n_regs = integer(),        # total no. of reg. models (logistic regressions)
    parfit_allowed = FALSE,    # allow parallel fit of multivar outvar when 1) reg$parfit = TRUE & 2) all.outvar.bin = TRUE
    initialize = function(reg, no_set_outvar = FALSE, ...) {
      self$reg <- reg
      if (!no_set_outvar) self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      # Number of sep. regressions to run, based on the number of outcomes in SingleRegressionFormClass length(reg$outvar)
      # or the number of objects in ListOfRegressionForms
      self$n_regs <- get_n_regs(reg)
      # self$n_regs <- length(reg$outvar) # switched to S3 dispatch
      all.outvar.bin <-  all(reg$outvar.class %in% gvars$sVartypes$bin)

      if (reg$parfit & all.outvar.bin & (self$n_regs > 1)) self$parfit_allowed <- TRUE
      #**** NOTE: for ltmle this should be changed to: if (reg$sep_predvars_sets) self$parfit_allowed <- TRUE
      if (gvars$verbose) {
        print("#----------------------------------------------------------------------------------")
        print("New instance of SummariesModel:")
        print("#----------------------------------------------------------------------------------")
        # if (reg$sep_predvars_sets) {
        if ("ListOfRegressionForms" %in% class(reg$RegressionForms)) {
          print("...ListOfRegressionForms..."); if (!is.null(names(reg$RegressionForms))) {print(names(reg$RegressionForms))}
          # print("Outcomes: "); print(unlist(lapply(reg$outvar, function(outvars) "(" %+% paste(outvars, collapse = ", ") %+% ")")))
          # print("Predictors: "); print(unlist(lapply(reg$predvars, function(predvars) "(" %+% paste(predvars, collapse = ", ") %+% ")")))
        } else {
          # reg$show()
          print("Outcomes: " %+% paste(reg$RegressionForms$outvar, collapse = ", "))
          print("Predictors: " %+% paste(reg$RegressionForms$predvars, collapse = ", "))
        }
        print("No. of regressions: " %+% self$n_regs)
        print("All outcomes binary? " %+% all.outvar.bin)
        # print("reg$subset_vars"); print(reg$subset_vars)
        # print("reg$subset_exprs"); print(reg$subset_exprs)
        # print("reg$outvar.class"); print(reg$outvar.class)
        if (self$parfit_allowed) print("Performing parallel fits: " %+% self$parfit_allowed)
        print("#----------------------------------------------------------------------------------")
      }

      # Factorize the joint into univariate regressions, by dimensionality of the outcome variable (sA_nms):
      for (k_i in 1:self$n_regs) {
        reg_i <- reg$clone()
        reg_i$ChangeManyToOneRegresssion(k_i, reg)
        # Calling the constructor for the summary model P(sA[j]|\bar{sA}[j-1], sW}), dispatching on reg_i class
        PsAsW.model <- newsummarymodel(reg = reg_i, ...)
        private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
        names(private$PsAsW.models)[k_i] <- "P(sA|sW)."%+%k_i
      }
      invisible(self)
    },
    length = function(){ base::length(private$PsAsW.models) },
    getPsAsW.models = function() { private$PsAsW.models },  # get all summary model objects (one model object per outcome var sA[j])
    getcumprodAeqa = function() { private$cumprodAeqa },  # get joint prob as a vector of the cumulative prod over j for P(sA[j]=a[j]|sW)
    fit = function(data) {
      assert_that(is.DataStorageClass(data))
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$fit(data = data)
        }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each fitRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        fitRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$fit(data = data)
        }
        # copy the fits one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.fit(fitRes[[k_i]])
        }
      }
      invisible(self)
    },
    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DataStorageClass(newdata))
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        }
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each predRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        predRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predict(newdata = newdata)
        }
        # copy the predictions one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.predict(predRes[[k_i]])
        }
      }
      invisible(self)
    },
    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Uses daughter objects (stored from prev call to fit()) to get predictions for P(sA=obsdat.sA|sW=sw)
    # Invisibly returns the joint probability P(sA=sa|sW=sw), also saves it as a private field "cumprodAeqa"
    # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    predictAeqa = function(newdata, ...) {
      assert_that(!missing(newdata))
      assert_that(is.DataStorageClass(newdata))
      n <- newdata$nobs
      if (!self$parfit_allowed) {
        cumprodAeqa <- rep.int(1, n)
        # loop over all regressions in PsAsW.models:
        for (k_i in seq_along(private$PsAsW.models)) {
          cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, ...)
        }
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = TRUE)
        probAeqa_list <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, ...)
        }
        probAeqa_mat <- do.call('cbind', probAeqa_list)
        cumprodAeqa <- matrixStats::rowProds(probAeqa_mat)
      }
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },

    sampleA = function(newdata, ...) {
      assert_that(!missing(newdata))
      assert_that(is.DataStorageClass(newdata))
      n <- newdata$nobs
      # loop over all regressions in PsAsW.models, sample CONDITIONALLY on observations that haven't been put in a specific bin yet
      sampleA_mat <- matrix(0L, nrow = n, ncol = length(private$PsAsW.models))
      for (k_i in seq_along(private$PsAsW.models)) {
        sampleA_newcat <- private$PsAsW.models[[k_i]]$sampleA(newdata = newdata, ...)
        if (k_i == 1L) sampleA_mat[, k_i] <- sampleA_newcat
        # carry forward all previously sampled 1's (degenerate ones a bin a chosen for the first time)
        if (k_i > 1) {
          # if you succeeded at the previous bin, your 1L is carried through till the end:
          sampleA_mat[(sampleA_mat[, k_i - 1] == 1L), k_i] <- 1L
          # if you haven't succeeded at the previous bin, you get a chance to succeed at this category:
          sampleA_mat[(sampleA_mat[, k_i - 1] == 0L), k_i] <- sampleA_newcat[(sampleA_mat[, k_i - 1] == 0L)]
        }
      }
      if (length(private$PsAsW.models) > 1) {
        sampleA_mat[, length(private$PsAsW.models)] <- 1L # make last category a reference category
        sampleA_cat <- rowSums(1L - sampleA_mat) + 1L
      } else {
        sampleA_cat <- as.vector(sampleA_mat)
      }
      return(sampleA_cat)
    },

    # call itself until reaches a terminal model fit with coefficients + regression returned with show()
    get.fits = function(format_table = FALSE) {
      res_models <- NULL
      for (k_i in seq_along(private$PsAsW.models)) {
        res <- private$PsAsW.models[[k_i]]$get.fits(format_table)
        if (is.list(res)) res_models <- c(res_models, res)
      }
      return(res_models)
    }
  ),
  active = list(
    # recursively call all saved daughter model fits and wipe out any traces of saved data
    wipe.alldat = function() {
      for (k_i in seq_along(private$PsAsW.models)) {
        private$PsAsW.models[[k_i]]$wipe.alldat
      }
      return(self)
    }
  ),
  private = list(
    deep_clone = function(name, value) {
      # if value is is an environment, quick way to copy:
      # list2env(as.list.environment(value, all.names = TRUE), parent = emptyenv())
      # if a list of R6 objects, make a deep copy of each:
      if (name == "PsAsW.models") {
        lapply(value, function(PsAsW.model) PsAsW.model$clone(deep=TRUE))
      # to check the value is an R6 object:
      } else if (inherits(value, "R6")) {
        value$clone(deep=TRUE)
      } else {
        # For all other fields, just return the value
        value
      }
    },
    PsAsW.models = list(),
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)


# -------------------------------------------------------------------------------------------
# Same code called from ContinSummaryModel$new and CategorSummaryModel$new:
# From a single categorical/continous outcome regression define a regression class for many binary dummy outcomes.
# First makes a clone of the parent RegressionClass and the recents the previous outvar/outvar.class to new binary outcomes/classes.
# Defines subset_var evaluation for new bins (for fitting the hazard of each new category/dummy)
# NOTE: This subsetting is performed by var name only (which automatically evaluates as !gvars$misval(var)) for speed & memory efficiency
# The code in DataStorageClass$binirize.sVar() will automatically set the indicator Bin_K[i] to NA when Bin_K-1[i] is 1 for the first time.
# Thus, by excluding all observations such that !is.na(Bin_K) we end up fitting the hazard for Bin_K=1.
# -------------------------------------------------------------------------------------------
def_regs_subset <- function(self) {
  bin_regs <- self$reg$clone()  # Instead of defining new RegressionClass, just clone the parent reg object and adjust the outcomes
  bin_regs$reg_hazard <- TRUE   # Don`t add degenerate bins as predictors in each binary regression
  if (!self$reg$pool_cont) {
    add.oldsubset <- TRUE
    new.subsets <- lapply(self$reg$bin_nms,
                              function(var) {
                                res <- var
                                if (add.oldsubset) res <- c(res, self$reg$subset_vars)
                                res
                              })

    new.sAclass <- as.list(rep_len(gvars$sVartypes$bin, self$reg$nbins))
    names(new.sAclass) <- self$reg$bin_nms
    bin_regs$outvar.class <- new.sAclass
    bin_regs$outvar <- self$reg$bin_nms
    bin_regs$predvars <- self$reg$predvars
    bin_regs$subset_vars <- new.subsets
    # bin_regs$ChangeOneToManyRegresssions(regs_list = list(outvar.class = new.sAclass,
    #                                                       outvar = self$reg$bin_nms,
    #                                                       predvars = self$reg$predvars,
    #                                                       subset_vars = new.subsets
    #                                                       ))
  # Same but when pooling across bin indicators:
  } else {
    bin_regs$outvar.class <- gvars$sVartypes$bin
    bin_regs$outvar <- self$outvar
    bin_regs$outvars_to_pool <- self$reg$bin_nms
    if (gvars$verbose)  {
      print("pooled bin_regs$outvar: "); print(bin_regs$outvar)
      print("bin_regs$outvars_to_pool: "); print(bin_regs$outvars_to_pool)
      print("bin_regs$subset_vars: "); print(bin_regs$subset_vars)
    }
  }
  bin_regs$resetS3class()
  return(bin_regs)
}

# -------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------
#' R6 class for fitting and predicting joint probability for a univariate continuous summary measure sA[j]
#'
#' This R6 class defines and fits a conditional probability model \code{P(sA[j]|sW,...)} for a univariate
#'  continuous summary measure \code{sA[j]}. This class inherits from \code{\link{SummariesModel}} class.
#'  Defines the fitting algorithm for a regression model \code{sA[j] ~ sW + ...}.
#'  Reconstructs the likelihood \code{P(sA[j]=sa[j]|sW,...)} afterwards.
#'  Continuous \code{sA[j]} is discretized using either of the 3 interval cutoff methods,
#'  defined via \code{\link{RegressionClass}} object \code{reg} passed to this class constructor.
#'  The fitting algorithm estimates the binary regressions for hazard \code{Bin_sA[j][i] ~ sW},
#'  i.e., the probability that continuous \code{sA[j]} falls into bin \code{i}, \code{Bin_sA[j]_i},
#'  given that \code{sA[j]} does not belong to any prior bins \code{Bin_sA[j]_1, ..., Bin_sA[j]_{i-1}}.
#'  The dataset of discretized summary measures (\code{BinsA[j]_1,...,BinsA[j]_M}) is created
#'  inside the passed \code{data} or \code{newdata} object while discretizing \code{sA[j]} into \code{M} bins.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{intrvls}} - .
#' \item{\code{intrvls.width}} - .
#' \item{\code{bin_weights}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, DataStorageClass.g0, DataStorageClass.gstar, ...)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
ContinSummaryModel <- R6Class(classname = "ContinSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the continous outcome var (sA[j])
    intrvls = NULL,
    intrvls.width = NULL,
    bin_weights = NULL,
    # Define settings for fitting contin sA and then call $new for super class (SummariesModel)
    initialize = function(reg, DataStorageClass.g0, DataStorageClass.gstar, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      if (is.null(reg$intrvls)) {
        assert_that(is.DataStorageClass(DataStorageClass.g0))
        self$intrvls <- DataStorageClass.g0$detect.sVar.intrvls(reg$outvar,
                                                      nbins = self$reg$nbins,
                                                      bin_bymass = self$reg$bin_bymass,
                                                      bin_bydhist = self$reg$bin_bydhist,
                                                      max_nperbin = self$reg$max_nperbin)
        if (!missing(DataStorageClass.gstar)) {
          assert_that(is.DataStorageClass(DataStorageClass.gstar))
          gstar.intrvls <- DataStorageClass.gstar$detect.sVar.intrvls(reg$outvar,
                                                      nbins = self$reg$nbins,
                                                      bin_bymass = self$reg$bin_bymass,
                                                      bin_bydhist = self$reg$bin_bydhist,
                                                      max_nperbin = self$reg$max_nperbin)
          self$intrvls <- unique(sort(union(self$intrvls, gstar.intrvls)))
        }
        # Define the number of bins (no. of binary regressions to run),
        # new outvar var names (bin names); all predvars remain unchanged;
        self$reg$intrvls <- self$intrvls
      } else {
        self$intrvls <- self$reg$intrvls
      }
      self$reg$nbins <- length(self$intrvls) - 1
      self$reg$bin_nms <- DataStorageClass.g0$bin.nms.sVar(reg$outvar, self$reg$nbins)
      # Save bin widths in reg class (naming the vector entries by bin names):
      self$intrvls.width <- diff(self$intrvls)
      self$intrvls.width[self$intrvls.width <= gvars$tolerr] <- 1
      self$reg$intrvls.width <- self$intrvls.width
      names(self$reg$intrvls.width) <- names(self$intrvls.width) <- self$reg$bin_nms
      if (gvars$verbose)  {
        print("ContinSummaryModel outcome: "%+%self$outvar)
      }
      bin_regs <- def_regs_subset(self = self)
      super$initialize(reg = bin_regs, no_set_outvar = TRUE, ...)
    },
    # Transforms data for continous outcome to discretized bins sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      assert_that(is.DataStorageClass(data))
      # Binirizes & saves binned matrix inside DataStorageClass
      data$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      if (gvars$verbose) {
        print("performing fitting for continuous outcome: " %+% self$outvar)
        print("freq counts by bin for continuous outcome: "); print(table(data$ord.sVar))
        print("binned dataset: "); print(head(cbind(data$ord.sVar, data$dat.bin.sVar), 5))
      }
      super$fit(data) # call the parent class fit method
      if (gvars$verbose) message("fit for outcome " %+% self$outvar %+% " succeeded...")
      data$emptydat.bin.sVar # wiping out binirized mat in data DataStorageClass object...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      invisible(self)
    },
    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DataStorageClass(newdata))
      if (gvars$verbose) print("performing prediction for continuous outcome: " %+% self$outvar)
      # mat_bin doesn't need to be saved (even though its invisibly returned); mat_bin is automatically saved in datnet.sW.sA - a potentially dangerous side-effect!!!
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      super$predict(newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DataStorageClass object...
      invisible(self)
    },
    # Convert contin. sA vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    predictAeqa = function(newdata) { # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a`s)
      assert_that(is.DataStorageClass(newdata))
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$reg$nbins, bin.nms = self$reg$bin_nms)
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      bws <- newdata$get.sVar.bw(name.sVar = self$outvar, intervals = self$intrvls)
      self$bin_weights <- (1 / bws) # weight based on 1 / (sVar bin widths)
      # Option 1: ADJUST FINAL PROB by bw.j TO OBTAIN density at a point f(sa|sw) = P(sA=sa|sW=sw):
      cumprodAeqa <- super$predictAeqa(newdata = newdata) * self$bin_weights
      # Alternative 2: ALso integrate the difference of sA value and its left most bin cutoff: x - b_{j-1} and pass it
      # This is done so that we can integrate the constant hazard all the way to the value of x:
        # * (1 - bw.j.sA_diff*(1/self$bin_weights)*probA1) (discrete)
        # * exp(-bw.j.sA_diff*(1/self$bin_weights)*probA1) (continuous)
      # bw.j.sA_diff <- newdata$get.sVar.bwdiff(name.sVar = self$outvar, intervals = self$intrvls)
      # cumprodAeqa <- super$predictAeqa(newdata = newdata, bw.j.sA_diff = bw.j.sA_diff) * self$bin_weights
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      self$bin_weights <- NULL # wiping out self$bin_weights...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },
    sampleA = function() {
      stop("not implemented")
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)

## ---------------------------------------------------------------------
#' R6 class for fitting and predicting joint probability for a univariate categorical summary measure sA[j]
#'
#' This R6 class defines and fits a conditional probability model \code{P(sA[j]|sW,...)} for a univariate
#'  categorical summary measure \code{sA[j]}. This class inherits from \code{\link{SummariesModel}} class.
#'  Defines the fitting algorithm for a regression model \code{sA[j] ~ sW + ...}.
#'  Reconstructs the likelihood \code{P(sA[j]=sa[j]|sW,...)} afterwards.
#'  Categorical \code{sA[j]} is first redefined into \code{length(levels)} bin indicator variables, where
#'  \code{levels} is a numeric vector of all unique categories in \code{sA[j]}.
#'  The fitting algorithm estimates the binary regressions for hazard for each bin indicator, \code{Bin_sA[j][i] ~ sW},
#'  i.e., the probability that categorical \code{sA[j]} falls into bin \code{i}, \code{Bin_sA[j]_i},
#'  given that \code{sA[j]} does not fall in any prior bins \code{Bin_sA[j]_1, ..., Bin_sA[j]_{i-1}}.
#'  The dataset of bin indicators (\code{BinsA[j]_1,...,BinsA[j]_M}) is created
#'  inside the passed \code{data} or \code{newdata} object when defining \code{length(levels)} bins for \code{sA[j]}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{levels}} - .
#' \item{\code{nbins}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, DataStorageClass.g0, ...)}}{...}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
CategorSummaryModel <- R6Class(classname = "CategorSummaryModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the categorical outcome var (sA[j])
    levels = numeric(),       # all unique values for sA[j] sorted in increasing order
    nbins = integer(),
    # Define settings for fitting cat sA and then call $new for super class (SummariesModel)
    initialize = function(reg, DataStorageClass.g0, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      # Define the number of bins (no. of binary regressions to run) based on number of unique levels for categorical sVar:
      # All predvars remain unchanged
      # browser()

      if (self$reg$get.reg$censoring & gvars$verbose) {
        message("...fitting a model for categorical censoring...")
      }

      if (is.null(reg$levels)) {
        assert_that(is.DataStorageClass(DataStorageClass.g0))
        self$levels <- self$reg$levels <- DataStorageClass.g0$detect.cat.sVar.levels(reg$outvar)
      } else {
        self$levels <- self$reg$levels
      }
      self$nbins <- self$reg$nbins <- length(self$levels)
      self$reg$bin_nms <- DataStorageClass.g0$bin.nms.sVar(reg$outvar, self$reg$nbins)
      if (gvars$verbose)  {
        print("CategorSummaryModel outcome: "%+%self$outvar)
      }
      bin_regs <- def_regs_subset(self = self)
      super$initialize(reg = bin_regs, no_set_outvar = TRUE, ...)
    },

    # Transforms data for categorical outcome to bin indicators sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      assert_that(is.DataStorageClass(data))
      # Binirizes & saves binned matrix inside DataStorageClass for categorical sVar
      data$binirize.sVar(name.sVar = self$outvar, levels = self$levels)
      if (gvars$verbose) {
        print("performing fitting for categorical outcome: " %+% self$outvar)
        print("freq counts by bin for categorical outcome: "); print(table(data$get.sVar(self$outvar)))
        print("binned dataset: "); print(head(cbind(sA = data$get.sVar(self$outvar), data$dat.bin.sVar), 5))
      }
      super$fit(data) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")
      data$emptydat.bin.sVar # wiping out binirized mat in data DataStorageClass object...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      invisible(self)
    },
    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DataStorageClass(newdata))
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      newdata$binirize.sVar(name.sVar = self$outvar, levels = self$levels)
      super$predict(newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DataStorageClass object...
      invisible(self)
    },
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    predictAeqa = function(newdata) {
      assert_that(is.DataStorageClass(newdata))
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      newdata$binirize.sVar(name.sVar = self$outvar, levels = self$levels)
      cumprodAeqa <- super$predictAeqa(newdata = newdata)
      newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      self$wipe.alldat # wiping out all data traces in ContinSummaryModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },
    sampleA = function(newdata) {
      assert_that(is.DataStorageClass(newdata))
      sampleA <- self$levels[super$sampleA(newdata = newdata)] # bring the sampled variable back to its original scale / levels:
      return(sampleA)
    }
  ),
  active = list(
    cats = function() {seq_len(self$reg$nbins)}
  )
)

StratifySummariesModel <- R6Class(classname = "StratifySummariesModel",
  inherit = SummariesModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the categorical outcome var (sA[j])
    # levels = numeric(),       # all unique values for sA[j] sorted in increasing order
    nbins = integer(),
    subset_exprs = NULL,
    # Define settings for fitting cat sA and then call $new for super class (SummariesModel)
    # We produce a regression class with 3 outvars (same) and 3 outvar.class (same)
    # The daughter regression clases resulting from this need to be of class StratifiedRegressionModelClass
    initialize = function(reg, DataStorageClass.g0, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      self$subset_exprs <- reg$subset_exprs
      assert_that(length(self$reg$subset_exprs) > 1L)
      assert_that(length(self$reg$outvar) == 1L)
      if (gvars$verbose)  {
        print("StratifySummariesModel outcome: "%+%self$outvar)
        print("StratifySummariesModel expressions: ("%+% paste(self$subset_exprs, collapse=",") %+% ")")
      }
      # print("self$reg before stratification:"); self$reg$show()
      stratify_regs <- self$reg$clone()  # Instead of defining new RegressionClass, just clone the parent reg object and adjust the outcomes
      stratify_regs$outvar <- rep_len(self$reg$outvar, length(self$reg$subset_exprs))
      stratify_regs$outvar.class <- as.list(rep_len(self$reg$outvar.class, length(self$reg$subset_exprs)))
      names(stratify_regs$outvar.class) <- stratify_regs$outvar
      # by turning stratify_regs$subset_exprs into a list the subsetting will be performed on subset_exprs during next S3 dispatch on SummaryModel
      stratify_regs$subset_exprs <- as.list(stratify_regs$subset_exprs)
      names(stratify_regs$subset_exprs) <- stratify_regs$outvar
      stratify_regs$reg_hazard <- TRUE
      # print("stratify_regs class:"); stratify_regs$show()
      super$initialize(reg = stratify_regs, no_set_outvar = TRUE, DataStorageClass.g0 = DataStorageClass.g0, ...)
    },

    # Transforms data for categorical outcome to bin indicators sA[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in sA - names have changed though)
    fit = function(data) {
      assert_that(is.DataStorageClass(data))
      if (gvars$verbose) {
        print("performing fitting for outcome based on stratified model for outcome: " %+% self$outvar)
        # print("following subsets are defined: "); print(table(data$get.sVar(self$outvar)))
      }
      super$fit(data) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")
      invisible(self)
    },
    # P(A^s=1|W^s=w^s): uses private$m.fit to generate predictions
    predict = function(newdata) {
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DataStorageClass(newdata))
      if (gvars$verbose) print("performing prediction for outcome based on stratified model: " %+% self$outvar)
      super$predict(newdata)
      invisible(self)
    },
    # Invisibly return cumm. prob P(sA=sa|sW=sw)
    # P(A^s=a^s|W^s=w^s) - calculating the likelihood for obsdat.sA[i] (n vector of a's):
    predictAeqa = function(newdata) {
      assert_that(is.DataStorageClass(newdata))
      if (gvars$verbose) print("performing prediction for outcome based on stratified model: " %+% self$outvar)
      cumprodAeqa <- super$predictAeqa(newdata = newdata)
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    },
    sampleA = function(newdata) {
      assert_that(is.DataStorageClass(newdata))
      # bring the sampled variable back to its original scale / levels:
      sampleA <- super$sampleA(newdata = newdata)
      return(sampleA)
    }
  ),
  active = list(
    # cats = function() {seq_len(self$reg$nbins)}
  )
)
