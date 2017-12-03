
# ---------------------------------------------------------------------------------------
#' Fit sequential (iterative means) GCOMP by pooling across more than one regimen
#'
#' Calls \code{fit_GCOMP} with argument \code{TMLE = FALSE}.
#' @param OData Input data object created by \code{importData} function.
#' @param intervened_TRT Column name in the input data with the probabilities (or indicators) of counterfactual treatment nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on treatment node(s).
#' @param ... Arguments that will be passed down to the underlying function \code{fit_GCOMP}
#' @return \code{data.table} with TMLE survival by time
#' @seealso \code{\link{fit_GCOMP}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fit_pooled_GCOMP <- function(OData, intervened_TRT, ...) {
  if (length(intervened_TRT) == 1L) stop("pooled GCOMP requires more than one intervention, for a single intervention please call 'fit_GCOMP' directly.")

  nodes <- OData$nodes
  copied_dat <- list()
  for (TRT in intervened_TRT) {
    Odat_DT_i <- copy(OData$dat.sVar[])
    Odat_DT_i[, "rule.name" := TRT][, "A.star" := eval(as.name(TRT))][, nodes$IDnode := paste0(eval(as.name(nodes$IDnode)), rule.name)]
    copied_dat[[TRT]] <- Odat_DT_i
  }

  Odat_DT_pooled <- rbindlist(copied_dat)

  ## todo: need a smarter way to just clone then copy OData R6 object, rather than having to re-import all the variables
  OData_pooled <- stremr::importData(Odat_DT_pooled,
                              ID = nodes$IDnode,
                              t = nodes$tnode,
                              covars = nodes$Lnodes,
                              CENS = nodes$Cnodes,
                              TRT = nodes$Anodes,
                              MONITOR = nodes$Nnodes,
                              OUTCOME = nodes$Ynode)

  # fit_GCOMP(TMLE = FALSE, ...)
  gcomp <- fit_GCOMP(intervened_TRT = c("A.star"),
                     OData = OData_pooled,
                     ...)
  return(gcomp)
}

# ---------------------------------------------------------------------------------------
#' Fit TMLE by pooling across more than one regimen
#'
#' Calls \code{fit_GCOMP} with argument \code{TMLE = FALSE}.
#' @param OData Input data object created by \code{importData} function.
#' @param intervened_TRT Column name in the input data with the probabilities (or indicators) of counterfactual treatment nodes being equal to 1 at each time point.
#' Leave the argument unspecified (\code{NULL}) when not intervening on treatment node(s).
#' @param IPWeights (Optional) result of calling function \code{getIPWeights} for running TMLE (evaluated automatically when missing)
#' @param ... Arguments that will be passed down to the underlying function \code{fit_GCOMP}
#' @return \code{data.table} with TMLE survival by time
#' @seealso \code{\link{fit_GCOMP}}
#' @example tests/examples/2_building_blocks_example.R
#' @export
fit_pooled_TMLE <- function(OData, intervened_TRT, IPWeights = NULL, ...) {
  if (length(intervened_TRT) == 1L) stop("pooled TMLE requires more than one intervention, for a single intervention please call 'fit_TMLE' directly.")
  nodes <- OData$nodes
  ## ----------------------------------------------------------------
  ## Evaluate the weights for each regimen, separately. Then rbind into a single dataset of weights.
  ## ----------------------------------------------------------------
  if (is.null(IPWeights)) {
    IPWeights <- purrr::map(intervened_TRT, getIPWeights, OData = OData)
  }
  IPWeights <- purrr::map(IPWeights, ~ .x[, nodes$IDnode := paste0(eval(as.name(nodes$IDnode)), rule.name)])
  IPWeights <- rbindlist(IPWeights)
  setkeyv(IPWeights, c(nodes$IDnode, nodes$tnode))
  attributes(IPWeights)[['intervened_TRT']] <- "A.star"
  attributes(IPWeights)[['intervened_MONITOR']] <- NULL

  fit_pooled_GCOMP(OData = OData, intervened_TRT = intervened_TRT, IPWeights = IPWeights, TMLE = TRUE, ...)
}
