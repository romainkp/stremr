# nocov start
#' Plot survival estimates using base R graphics
#' @param surv_list A list with survival estimates, one per regimen.
#' @param t The vector of time values for plotting.
#' @param t_int_sel Optional vector of indices that subsets \code{t}.
#' If omitted the survival for all t values will be plotted.
#' @param y_lab y-axis title.
#' @param x_lab x-axis title.
#' @param miny Minimum y value to plot
#' @param x_legend y-coordinate for legend location.
#' @param y_legend x-coordinate for legend location.
#' @param cex Same as R plot function.
#' @param ... Additional arguments to be passed on to base R plot function.
#' @export
#'
f_plot_survest <- function(surv_list, t, t_int_sel, y_lab, x_lab, miny, x_legend, y_legend, cex = 0.7, ...) {
  # ptsize <- 0.7
  # ptsize <- 0.4
  warning("Deprecated: f_plot_survest is now considered deprecated and will be removed from stremr in the forthcoming releases.
Consider using ggsurv and ggRD functions instead.")
  counter <- 0
  if (missing(y_lab)) y_lab <- ""
  if (missing(x_lab)) x_lab <- "Follow-up period since study entry"
  if (missing(t)) t <- seq_along(surv_list[[1]])
  if (missing(t_int_sel)) t_int_sel <- seq_along(t)
  if (missing(miny)) miny <- min(unlist(lapply(surv_list, function(x) min(x[t_int_sel], na.rm = TRUE))))
  # if (missing(x_legend)) x_legend <- (max(t_int_sel, na.rm = TRUE) - min(t_int_sel, na.rm = TRUE)) * 2/3 + min(t_int_sel, na.rm = TRUE)
  # if (missing(y_legend)) y_legend <- (1 - miny) * 4/5 + miny
  for(d.j in names(surv_list)){
    counter <- counter + 1
    plot(as.integer(t[t_int_sel]), surv_list[[d.j]][t_int_sel], col = counter, type = 'b', cex = cex, ylim = c(miny, 1), ylab = y_lab, xlab = x_lab)
    par(new=TRUE)
  }
  if (missing(y_legend)) {
    if (missing(x_legend)) {
      x_legend <- "bottomleft"
    } else {
      if (!is.character(x_legend)) stop("x_legend must be a character when y_legend is unspecified")
    }
    legend(x_legend, legend = names(surv_list), col = c(1:length(names(surv_list))), cex = cex, pch = 1)
  } else {
    legend(x_legend, y_legend, legend = names(surv_list), col = c(1:length(names(surv_list))), cex = cex, pch = 1)
  }
}

f_obtain_St <- function(sysArg, est_obj, optArgReport, est_name = "St.TMLE", t_name) {
  sysArg[["surv_list"]] <- lapply(est_obj, '[[', est_name)
  rule.names <- unlist(lapply(est_obj, function(est_obj_res) est_obj_res[['rule.name']][1]))
  names(sysArg[["surv_list"]]) <- rule.names
  sysArg[["t"]] <- est_obj[[1]][[t_name]]
  userArg <- intersect(names(formals(f_plot_survest)), names(optArgReport)) # captures optional arguments given by user for customizing report
  if(length(userArg) > 0) sysArg <- c(sysArg, optArgReport[userArg])
  return(sysArg)
}

# nocov end