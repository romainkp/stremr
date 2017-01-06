
#-----------------------------------------------------------------------------
# Global State Vars (can be controlled globally with options(stremr.optname = ))
#-----------------------------------------------------------------------------
gvars <- new.env(parent = emptyenv())
gvars$verbose <- FALSE      # verbose mode (print all messages)
gvars$opts <- list()        # named list of package options that is controllable by the user (set_all_stremr_options())
gvars$misval <- NA_integer_ # the default missing value for observations (# gvars$misval <- -.Machine$integer.max)
gvars$misXreplace <- 0L     # the default replacement value for misval that appear in the design matrix
gvars$tolerr <- 10^-12      # tolerance error: assume for abs(a-b) < gvars$tolerr => a = b
gvars$sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
gvars$noCENScat <- 0L       # the reference category that designates continuation of follow-up

allowed.fit.package <- c("speedglm", "glm", "h2o")
allowed.fit.algorithm = c("glm", "gbm", "randomForest", "deeplearning")
allowed.bin.method = c("equal.mass", "equal.len", "dhist")

#' Querying/setting a single \code{stremr} option
#'
#' To list all \code{stremr} options, just run this function without any parameters provided. To query only one value, pass the first parameter. To set that, use the \code{value} parameter too.
#'
#' The arguments of \code{\link{set_all_stremr_options}} list all available \code{stremr} options.
#'
#' @param o Option name (string). See \code{\link{set_all_stremr_options}}.
#' @param value Value to assign (optional)
#' @export
#' @seealso \code{\link{set_all_stremr_options}}
#' @examples \dontrun{
#' stremrOptions()
#' stremrOptions('fit.package')
#' stremrOptions('fit.package', 'h2o')
#' }
stremrOptions <- function (o, value)  {
  res <- getOption("stremr")
  if (missing(value)) {
    if (missing(o))
        return(res)
    if (o %in% names(res))
        return(res[[o]])
    print("Possible `stremr` options:")
    print(names(res))
    stop(o %+% ": this options does not exist")
  } else {
    if (!o %in% names(res))
      stop(paste("Invalid option name:", o))
    if (is.null(value)) {
      res[o] <- list(NULL)
    }
    else {
      res[[o]] <- value
    }
    do.call("set_all_stremr_options", res)
  }
}

getopt <- function(optname) {
  return(stremrOptions(o = optname))
}

#' Print Current Option Settings for \code{stremr}
#' @return Invisibly returns a list of \code{stremr} options.
#' @seealso \code{\link{set_all_stremr_options}}
#' @export
print_stremr_opts <- function() {
  print(gvars$opts)
  invisible(gvars$opts)
}

#' Setting \code{stremr} Options
#'
#' Options that control \code{stremr} package.
#' \strong{Will reset all unspecified options (omitted arguments) to their default values}.
#' The preferred way to set options for \code{stremr} is to use \code{\link{stremrOptions}}, which allows specifying individual options without having to reset all other options.
#' To reset all options to their defaults simply run \code{set_all_stremr_options()} without any parameters/arguments.
#' @param fit.package Specify the default package for performing model fitting: c("speedglm", "glm", "h2o")
#' @param fit.algorithm Specify the default fitting algorithm: c("glm", "gbm", "randomForest", "deeplearning", "SuperLearner")
#' @param bin.method The method for choosing bins when discretizing and fitting the conditional continuous summary
#'  exposure variable \code{sA}. The default method is \code{"equal.len"}, which partitions the range of \code{sA}
#'  into equal length \code{nbins} intervals. Method \code{"equal.mass"} results in a data-adaptive selection of the bins
#'  based on equal mass (equal number of observations), i.e., each bin is defined so that it contains an approximately
#'  the same number of observations across all bins. The maximum number of observations in each bin is controlled
#'  by parameter \code{maxNperBin}. Method \code{"dhist"} uses a mix of the above two approaches,
#'  see Denby and Mallows "Variations on the Histogram" (2009) for more detail.
# @param parfit Default is \code{FALSE}. Set to \code{TRUE} to use \code{foreach} package and its functions
#  \code{foreach} and \code{dopar} to perform
#  parallel logistic regression fits and predictions for discretized continuous outcomes. This functionality
#  requires registering a parallel backend prior to running \code{stremr} function, e.g.,
#  using \code{doParallel} R package and running \code{registerDoParallel(cores = ncores)} for integer
#  \code{ncores} parallel jobs. For an example, see a test in "./tests/RUnit/RUnit_tests_04_netcont_sA_tests.R".
#' @param nbins Set the default number of bins when discretizing a continous outcome variable under setting
#'  \code{bin.method = "equal.len"}.
#'  If left as \code{NA} the total number of equal intervals (bins) is determined by the nearest integer of
#'  \code{nobs}/\code{maxNperBin}, where \code{nobs} is the total number of observations in the input data.
#' @param maxncats Max number of unique categories a categorical variable \code{sA[j]} can have.
#' If \code{sA[j]} has more it is automatically considered continuous.
# @param poolContinVar Set to \code{TRUE} for fitting a pooled regression which pools bin indicators across all bins.
# When fitting a model for binirized continuous outcome, set to \code{TRUE}
# for pooling bin indicators across several bins into one outcome regression?
#' @param maxNperBin Max number of observations per 1 bin for a continuous outcome (applies directly when
#'  \code{bin.method="equal.mass"} and indirectly when \code{bin.method="equal.len"}, but \code{nbins = NA}).
#' @param lower_bound_zero_Q Set to \code{TRUE} to bound the observation-specific Qs during the TMLE update step away from zero (with minimum value set at 10^-4).
#' Can help numerically stabilize the TMLE intercept estimates in some small-sample cases. Has no effect when \code{TMLE} = \code{FALSE}.
#' @param skip_update_zero_Q Set to \code{FALSE} to perform TMLE update with glm even when all of the Q's are zero.
#' When set to \code{TRUE} the TMLE update step is skipped if the predicted Q's are either all 0 or near 0, with TMLE intercept being set to 0.
#' @return Invisibly returns a list with old option settings.
#' @seealso \code{\link{stremrOptions}}, \code{\link{print_stremr_opts}}
#' @export
set_all_stremr_options <- function( fit.package = c("speedglm", "glm", "h2o"),
                            fit.algorithm = c("glm", "gbm", "randomForest", "deeplearning", "SuperLearner"),
                            bin.method = c("equal.mass", "equal.len", "dhist"),
                            nbins = 10,
                            maxncats = 20,
                            # poolContinVar = FALSE,
                            maxNperBin = 500,
                            lower_bound_zero_Q = TRUE,
                            skip_update_zero_Q = TRUE
                            ) {

  old.opts <- gvars$opts

  fit.package <- fit.package[1L]
  fit.algorithm <- fit.algorithm[1L]
  bin.method <- bin.method[1]
  if (!(fit.package %in% allowed.fit.package)) stop("fit.package must be one of: " %+% paste0(allowed.fit.package, collapse=", "))
  if (!(fit.algorithm %in% allowed.fit.algorithm)) stop("fit.algorithm must be one of: " %+% paste0(allowed.fit.algorithm, collapse=", "))
  if (!(bin.method %in% allowed.bin.method)) stop("bin.method must be one of: " %+% paste0(allowed.bin.method, collapse=", "))

  opts <- list(
    fit.package = fit.package,
    fit.algorithm = fit.algorithm,
    bin.method = bin.method,
    # parfit = parfit,
    nbins = nbins,
    maxncats = maxncats,
    # poolContinVar = poolContinVar,
    maxNperBin = maxNperBin,
    lower_bound_zero_Q = lower_bound_zero_Q,
    skip_update_zero_Q = skip_update_zero_Q
  )
  gvars$opts <- opts
  options(stremr = opts)
  invisible(old.opts)
}

# returns a function (alternatively a call) that tests for missing values in (sA, sW)
testmisfun <- function() {
  if (is.na(gvars$misval)) {
    return(is.na)
  } else if (is.null(gvars$misval)){
    return(is.null)
  } else if (is.integer(gvars$misval)) {
    return(function(x) {x==gvars$misval})
  } else {
    return(function(x) {x%in%gvars$misval})
  }
}

get.misval <- function() {
  gvars$misfun <- testmisfun()
  gvars$misval
}

set.misval <- function(gvars, newmisval) {
  oldmisval <- gvars$misval
  gvars$misval <- newmisval
  gvars$misfun <- testmisfun()    # EVERYTIME gvars$misval HAS CHANGED THIS NEEDS TO BE RESET/RERUN.
  invisible(oldmisval)
}
gvars$misfun <- testmisfun()

# Allows stremr functions to use e.g., getOption("stremr.verbose") to get verbose printing status
.onLoad <- function(libname, pkgname) {
  # reset all options to their defaults on load:
  set_all_stremr_options()
  op <- options()
  op.stremr <- list(
    stremr.verbose = gvars$verbose,
    stremr.file.path = tempdir(),
    # stremr.file.name = 'stremr-report-%T-%N-%n'
    stremr.file.name = 'stremr-report-'%+%Sys.Date()
  )
  toset <- !(names(op.stremr) %in% names(op))
  if (any(toset)) options(op.stremr[toset])
  invisible()
}

# Runs when attached to search() path such as by library() or require()
.onAttach <- function(...) {
  if (interactive()) {
  	packageStartupMessage('stremr')
  	# packageStartupMessage('Version: ', utils::packageDescription('stremr')$Version)
  	packageStartupMessage('Version: ', utils::packageDescription('stremr')$Version, '\n')
  	packageStartupMessage(
  "stremr IS IN EARLY DEVELOPMENT STAGE.
Please be to sure to check for frequent updates and report bugs at: http://github.com/osofr/stremr
To install the latest development version of stremr, please type this in your terminal:
  devtools::install_github('osofr/stremr')", '\n')
  	# packageStartupMessage('To see the vignette use vignette("stremr_vignette", package="stremr"). To see all available package documentation use help(package = "stremr") and ?stremr.', '\n')
  	# packageStartupMessage('To see the latest updates for this version, use news(package = "stremr").', '\n')
  }
}