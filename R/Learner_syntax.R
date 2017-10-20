# @importFrom gridisl defModel
# @name defModel
# @rdname defModel
# @export
# NULL

#' S3 methods for printing a collection of learners
#'
#' Prints the stack models
#' @param x An object (list) of class ModelStack
#' @param ... Additional options passed on to \code{print.PredictionModel}.
#' @export
print.ModelStack <- function(x, ...) {
  str(x)
  return(invisible(NULL))
}

# ---------------------------------------------------------------------------------------
#' Interface for defining models
#'
#' Specify either a single model or a grid of multple models by invoking the optional argument \code{param_grid}.
#' @param estimator A character string name of package and estimator (algorithm) name, separated by "__".
#' @param x A vector containing the subset of the names of the predictor variables to use in building this
#' particular learner or grid. This argument can be used to over-ride the values of \code{x} provided to \code{fit} function.
#' As such, the names supplied here must always be a subset of the names specified to \code{fit}.
#' When this argument is missing (default) the column names provided to \code{fit} are used as predictors in building this model / grid.
#' @param search_criteria Search criteria
#' @param param_grid Named list of model hyper parameters (optional).
#' Each named item defines either a fixed model parameter value or a vector of possible values.
#' In the latter case this function call would define a grid (ensemble) of models.
#' Each model in the ensemble corresponds by a particular fixed parameter value.
#' @param ... Additional modeling parameters to be passed on directly to the modeling function.
#' @export
defModel <- function(estimator, x, search_criteria, param_grid, ...) {
  pkg_est <- strsplit(estimator, "__", fixed = TRUE)[[1]]
  pkg <- pkg_est[1]
  if (length(pkg_est) > 1) est <- pkg_est[2] else est <- NULL

  ## call outside fun that parses ... and checks all args are named
  sVar.exprs <- capture.exprs(...)

  GRIDparams = list(fit.package = pkg, fit.algorithm = "grid", grid.algorithm = est, estimator = estimator)
  if (!missing(x)) GRIDparams[["x"]] <- x
  if (!missing(search_criteria)) GRIDparams[["search_criteria"]] <- search_criteria
  if (!missing(param_grid)) {
    if (!is.list(param_grid)) stop("'param_grid' must be a named list of model hyper parameters")
    if (is.null(names(param_grid)) || "" %in% names(param_grid)) stop("all items in 'param_grid' must be named")
    GRIDparams[["param_grid"]] <- param_grid
  }

  if (length(sVar.exprs) > 0) GRIDparams <- c(GRIDparams, sVar.exprs)
  GRIDparams <- list(GRIDparams)
  class(GRIDparams) <- c(class(GRIDparams), "ModelStack")
  return(GRIDparams)
}

# S3 method '+' for adding two ModelStack objects
# Summary measure lists in both get added as c(,) into the summary measures in model1 / model2 objects
#' @rdname defModel
#' @param model1 An object returned by a call to \code{defModel} function.
#' @param model2 An object returned by a call to \code{defModel} function.
#' @export
`+.ModelStack` <- function(model1, model2) {
  assert_that(is.ModelStack(model1))
  assert_that(is.ModelStack(model2))
  newStack <- append(model1, model2)
  class(newStack) <- c(class(newStack), "ModelStack")
  return(newStack)
}
