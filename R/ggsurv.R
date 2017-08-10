## -----------------------------------------------------------------------------
## Survival plots with ggplot2
## Functions adapted from GGally R package: https://github.com/ggobi/ggally
## Original R code by Edwin Thoen \email{edwinthoen@@gmail.com},
## modified by Oleg Sofrygin \email{oleg.sofrygin@@gmail.com}
## See http://ggobi.github.io/ggally/#ggallyggsurv for additional modifications of the resulting survival plot
## -----------------------------------------------------------------------------

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("cens", "surv", "CI95up", "CI95low", "up", "low", "time", "group",
                           "dx1", "dx2", "dx1_name", "dx2_name", "RD",
                           "RD.SE", "time", "contrast", "element_text", "estimates", "rank1", "rank2", "time_idx"))
}

#' Survival curves with ggplot2
#'
#' This function produces plots of the survival point estimates using \code{ggplot2}.
#' As a first argument it needs a list of estimates created by the one of the estimation function of the \code{stremr} package.
#' See http://ggobi.github.io/ggally/#ggallyggsurv for additional modifications of the resulting output plot.
#'
#' @export
#' @param estimates A list, one item per regime / intervention. Each list item must be a data.frame containing the
#' survival estimates by time for a single regime / intervention.
#' @param CI should a 95\% confidence interval be plotted? Defaults to \code{TRUE}.
#' Uses the standard error estimates provided as a separate column of the input data.
#' @param CI_line When \code{TRUE} the 95\% CIs will be plotted as a line function (same as main plot type).
#' When \code{FALSE} the 95\% CIs are plotted using \code{ggplot2::geom_ribbon}.
#' @param plot_cens mark the censored observations?
#' @param surv_col colour of the survival estimate. Defaults to black for
#'    one stratum, and to the default \code{ggplot2} colours for multiple
#'    strata. Length of vector with colour names should be either 1 or equal
#'    to the number of strata.
#' @param cens_col colour of the points that mark censored observations.
#' @param lty_est linetype of the survival curve(s). Vector length should be
#'    either 1 or equal to the number of strata.
#' @param shape_est shape type of the survival point estimates. Vector length should be
#'  either 1 or equal to the number of strata.
#' @param lty_ci linetype of the bounds that mark the 95\% CI.
#' @param size_est line width of the survival curve
#' @param size_ci line width of the 95\% CI
#' @param size_pt point size of the survival estimate at each time-point
#' @param cens_size point size of the censoring points
#' @param cens_shape shape of the points that mark censored observations.
#' @param back_white if \code{TRUE} the background will not be the default
#'    grey of \code{ggplot2}, but will be white with borders around the plot. Defaults to \code{TRUE}
#' @param xlab the label of the x-axis.
#' @param ylab the label of the y-axis.
#' @param main the plot label.
#' @param legend_pos Either the coordinates of the legend position inside the plot
#' (e.g., (0.9, 0.2)) or
#' the character word denoting the legend orientation with respect to the plot
#' (e.g., "bottom", "right" or "left").
#' @param surv_name The name of the column containing the survival estimates.
#' @param SE_name The name of the column containing the standard errors (SE) for each time-point estimate of survival.
#' @param order_legend Set to \code{TRUE} to order the legend display by final
#' survival time (highest first).
#' @param t_int_sel The subset of time-point indices for which survival should be plotted.
#' @param ymin The minimum value of the y axis. The default (\code{ymin=NULL}) is to use \code{ggplot}
#' to automatically adjust the limits of the y-axis.
#' @param ... Additional arguments (not used).
#' @return An object of class \code{ggplot}
#' @author Original R code by Edwin Thoen \email{edwinthoen@@gmail.com}, modified by Oleg Sofrygin \email{oleg.sofrygin@@gmail.com}
ggsurv <- function(
  estimates,
  CI         = TRUE,
  CI_line    = FALSE,
  plot_cens  = TRUE,
  surv_col   = 'gg.def',
  cens_col   = 'gg.def',
  lty_est    = 1,
  shape_est  = seq_along(estimates),
  lty_ci     = 2,
  size_est   = 0.5,
  size_ci    = size_est,
  size_pt    = size_est + 0.3,
  cens_size  = 2,
  cens_shape = 3,
  back_white = TRUE,
  xlab       = 'Time',
  ylab       = 'Survival',
  main       = '',
  legend_pos = "right",
  # surv_name = "St." %+% attr(estimates[[1]], "estimator_short"),
  # SE_name = "SE." %+% attr(estimates[[1]], "estimator_short"),
  surv_name = "St." %+% estimates[[1L]][["est_name"]][1L],
  SE_name = "SE." %+% estimates[[1L]][["est_name"]][1L],
  order_legend = TRUE,
  t_int_sel = NULL,
  ymin = NULL,
  ...
){

  if ("data.frame" %in% class(estimates)) {
    estimates <- estimates[[1]]
    shape_est <- seq_along(estimates)
    # surv_name <- "St." %+% attr(estimates[[1]], "estimator_short")
    # SE_name <- "SE." %+% attr(estimates[[1]], "estimator_short")
  }
  # if ("data.frame" %in% class(estimates)) estimates <- data.table::rbindlist(estimates[[1]])


  n.grps <- length(estimates)
  gr.name <- "regime"
  gr.df <- vector('list', n.grps)

  for (i in 1:n.grps) {
    surv_dat <- estimates[[i]]
    gr.df[[i]] <- data.table::data.table(
        time  = surv_dat[["time"]],
        surv  = surv_dat[[surv_name]],
        group = surv_dat[["rule.name"]]
      )
      if (SE_name %in% names(surv_dat)) {
        gr.df[[i]][, ("up") := surv + abs(qnorm(0.025))*surv_dat[[SE_name]]]
        gr.df[[i]][, ("low") := surv - abs(qnorm(0.025))*surv_dat[[SE_name]]]
      } else {
        CI <- FALSE
      }
      if (!is.null(t_int_sel)) gr.df[[i]] <- gr.df[[i]][t_int_sel, ]
  }

  dat      <- data.table::rbindlist(gr.df)

  pl <- ggplot2::ggplot(dat, ggplot2::aes(x = time, y = surv, group = group)) +
    ggplot2::geom_line(ggplot2::aes(col = group, lty = group), size = size_est) +
    ggplot2::geom_point(ggplot2::aes(col = group, shape = group), size = size_pt) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(main)

  pl <- if(surv_col[1] != 'gg.def'){
    scaleValues <- if (length(surv_col) == 1) {
      rep(surv_col, n.grps)
    } else{
      surv_col
    }
    pl +
      ggplot2::scale_colour_manual(name = gr.name, values = scaleValues) +
      ggplot2::scale_fill_manual(name = gr.name, values = scaleValues)
  } else {
    pl +
      ggplot2::scale_colour_discrete(name = gr.name) +
      ggplot2::scale_fill_discrete(name = gr.name)
  }

  lineScaleValues <- if (length(lty_est) == 1) {
    rep(lty_est, n.grps)
  } else {
    lty_est
  }
  pl <- pl + ggplot2::scale_linetype_manual(name = gr.name, values = lineScaleValues)

  pointShapeValues <- if (length(shape_est) == 1) {
    rep(shape_est, n.grps)
  } else {
    shape_est
  }
  pl <- pl + ggplot2::scale_shape_manual(name = gr.name, values = pointShapeValues)

  if(identical(CI,TRUE)) {
    if(length(surv_col) > 1 && length(lty_est) > 1){
      stop('Either surv_col or lty_est should be of length 1 in order to plot 95% CI with multiple n.grps')
    }

    stepLty <- if ((length(surv_col) > 1 | surv_col == 'gg.def')[1]) {
      lty_ci
    } else {
      surv_col
    }

    if (CI_line) {
      pl <- pl +
        ggplot2::geom_line(ggplot2::aes(y = up, lty = group, col = group), lty = stepLty, size = size_ci) +
        ggplot2::geom_line(ggplot2::aes(y = low,lty = group, col = group), lty = stepLty, size = size_ci)
    } else {
      pl <- pl +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = up, fill = group, linetype = group), alpha = 0.1, size = size_ci, lty = stepLty)
    }
  }

  if(identical(back_white, TRUE)) pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::theme(legend.position = legend_pos)
  if (!is.null(ymin)) pl <- pl + ggplot2::coord_cartesian(ylim = c(ymin, 1))
  return(pl)
}


#' Plot risk differences over time with ggplot2
#'
#' This function produces plots of risk differences over time using \code{ggplot2}.
#' As a first argument it needs a dataset of risk differences produced by \code{link{get_RDs}}.
#' See http://ggobi.github.io/ggally/#ggallyggsurv for additional modifications of the resulting output plot.
#'
#' @export
#' @param RDests A table in long format with a single column containing names of two contrasting regimens.
#' Each list item must be a data.frame containing the
#' risk difference estimates over time.
#' @param CI should a 95\% confidence interval be plotted? Defaults to \code{TRUE}.
#' Uses the standard error RDests provided as a separate column of the input data.
#' @param CI_errorbar Plot CIs as error bars around point estimates, default is \code{FALSE}.
#' @param CI_line When \code{TRUE} the 95\% CIs will be plotted as a line function (same as main plot type).
#' When \code{FALSE} the 95\% CIs are plotted using \code{ggplot2::geom_ribbon}.
#' @param plot_cens mark the censored observations?
#' @param surv_col colour of the survival estimate. Defaults to black for
#'    one stratum, and to the default \code{ggplot2} colours for multiple
#'    strata. Length of vector with colour names should be either 1 or equal
#'    to the number of strata.
#' @param cens_col colour of the points that mark censored observations.
#' @param lty_est linetype of the survival curve(s). Vector length should be
#'    either 1 or equal to the number of strata.
#' @param shape_est shape type of the survival point RDests. Vector length should be
#'  either 1 or equal to the number of strata.
#' @param lty_ci linetype of the bounds that mark the 95\% CI.
#' @param size_est line width of the survival curve
#' @param size_ci line width of the 95\% CI
#' @param size_pt point size of the survival estimate at each time-point
#' @param cens_size point size of the censoring points
#' @param cens_shape shape of the points that mark censored observations.
#' @param back_white if \code{TRUE} the background will not be the default
#'    grey of \code{ggplot2}, but will be white with borders around the plot. Defaults to \code{TRUE}
#' @param xlab the label of the x-axis.
#' @param ylab the label of the y-axis.
#' @param main the plot label.
#' @param legend_pos Either the coordinates of the legend position inside the plot
#' (e.g., (0.9, 0.2)) or
#' the character word denoting the legend orientation with respect to the plot
#' (e.g., "bottom", "right" or "left").
#' @param RD_name The name of the column containing the risk differences.
#' @param SE_name The name of the column containing the standard errors (SE) for each risk difference.
#' @param CIlow_name The name of the column containing the lower bound of the confidence interval (optional).
#' @param CIup_name The name of the column containing the upper bound of the confidence interval (optional).
#' @param order_legend Set to \code{TRUE} to order the legend display by final
#' survival time (highest first).
#' @param t_int_sel The subset of time-point indices for which survival should be plotted.
#' @param ymin The minimum value of the y axis. The default (\code{ymin=NULL}) is to use \code{ggplot}
#' to automatically adjust the limits of the y-axis.
#' @param ymax The maximum value of the y axis. The default (\code{ymax=NULL}) is to use \code{ggplot}
#' to automatically adjust the limits of the y-axis.
#' @param facet Set to \code{TRUE} to create a facet of plots by first / last tx name
#' @param font_size Font size for facet labels, legend text and legend title.
#' @param line_RD_plot Connect the point risk differences with a plotted line (using \code{geom_line}).
#' @param x_breaks Vector of breaks to display on x-axis
#' @param axis_font_size Control the font size of x and y axis, leave as NULL for default.
#' @param x_axis_font_angle Control the angle of rotation for x-axis labels (default is 0).
#' @param ... Additional arguments (not used).
#' @return An object of class \code{ggplot}
#' @author Original R code by Edwin Thoen \email{edwinthoen@@gmail.com}, modified by Oleg Sofrygin \email{oleg.sofrygin@@gmail.com}
ggRD <- function(
  RDests,
  CI         = TRUE,
  CI_errorbar = FALSE,
  CI_line    = FALSE,
  plot_cens  = TRUE,
  surv_col   = 'gg.def',
  cens_col   = 'gg.def',
  lty_est    = 1,
  shape_est  = seq_along(RDests),
  lty_ci     = 2,
  size_est   = 0.5,
  size_ci    = size_est,
  size_pt    = size_est + 0.3,
  cens_size  = 2,
  cens_shape = 3,
  back_white = TRUE,
  xlab       = 'Time',
  ylab       = 'Risk Difference',
  main       = '',
  legend_pos = "right",
  RD_name = "RD",
  SE_name = "RD.SE",
  CIlow_name = "CI95low",
  CIup_name = "CI95up",
  order_legend = TRUE,
  t_int_sel = NULL,
  ymin = NULL,
  ymax = NULL,
  facet = FALSE,
  font_size = 6,
  line_RD_plot = TRUE,
  x_breaks = unique(RDests[["time"]]),
  axis_font_size = NULL,
  x_axis_font_angle = 0,
  ...
){

  if (!"RD.SE" %in% names(RDests)) {
    RDests <- RDests %>% dplyr::mutate(RD.SE = NA)
  }

  dat <-  RDests %>%
          dplyr::rename_(RD = RD_name, RD.SE = SE_name)

  if (!is.null(t_int_sel))
    dat <- dat %>% dplyr::filter(time_idx %in% t_int_sel)

  n.grps <- length(unique(dat[["contrast"]]))
  gr.name <- "contrast"

  pl <- ggplot2::ggplot(dat, ggplot2::aes(x = time, y = RD, group = contrast))
  if (line_RD_plot) pl <- pl + ggplot2::geom_line(ggplot2::aes(col = contrast, lty = contrast), size = size_est)
  pl <- pl +
        # , shape = contrast
        ggplot2::geom_point(ggplot2::aes(col = contrast), size = size_pt) +
        # ggplot2::scale_x_discrete(breaks = waiver()) +
        ggplot2::scale_x_continuous(breaks = x_breaks) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::ggtitle(main)

  pl <- if(surv_col[1] != 'gg.def'){
    scaleValues <- if (length(surv_col) == 1) {
      rep(surv_col, n.grps)
    } else{
      surv_col
    }
    pl +
      ggplot2::scale_colour_manual(name = gr.name, values = scaleValues) +
      ggplot2::scale_fill_manual(name = gr.name, values = scaleValues)
  } else {
    pl +
      ggplot2::scale_colour_discrete(name = gr.name) +
      ggplot2::scale_fill_discrete(name = gr.name)
  }

  lineScaleValues <- if (length(lty_est) == 1) {
    rep(lty_est, n.grps)
  } else {
    lty_est
  }
  pl <- pl + ggplot2::scale_linetype_manual(name = gr.name, values = lineScaleValues)

  pointShapeValues <- if (length(shape_est) == 1) {
    rep(shape_est, n.grps)
  } else {
    shape_est
  }
  pl <- pl + ggplot2::scale_shape_manual(name = gr.name, values = pointShapeValues)

  if(identical(CI,TRUE)) {
    if(length(surv_col) > 1 && length(lty_est) > 1){
      stop('Either surv_col or lty_est should be of length 1 in order to plot 95% CI with multiple n.grps')
    }

    stepLty <- if ((length(surv_col) > 1 | surv_col == 'gg.def')[1]) {
      lty_ci
    } else {
      surv_col
    }

    if (CI_errorbar) {
      pl <- pl +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin = CIlow_name, ymax = CIup_name, lty = "contrast", col = "contrast"), lty = stepLty, size = size_ci)

    } else  if (CI_line) {
      pl <- pl +
        ggplot2::geom_line(ggplot2::aes_string(y = CIup_name, lty = "contrast", col = "contrast"), lty = stepLty, size = size_ci) +
        ggplot2::geom_line(ggplot2::aes_string(y = CIlow_name,lty = "contrast", col = "contrast"), lty = stepLty, size = size_ci)

    } else {
      pl <- pl +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = CIlow_name, ymax = CIup_name, fill = "contrast", linetype = "contrast"), alpha = 0.1, size = size_ci, lty = stepLty)
    }
  }

  if(identical(back_white, TRUE)) pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::theme(legend.position = legend_pos,
                            legend.text = element_text(size = font_size),
                            legend.title = element_text(size = font_size))
  if (!is.null(ymin) && !is.null(ymax)) pl <- pl + ggplot2::coord_cartesian(ylim = c(ymin, ymax))


  label_value2 <- function (labels, multi_line = TRUE) {
    labels <- lapply(labels, as.character)
    if (multi_line) {
      labels
    } else {
      out <- do.call("Map", c(list(paste, sep = " - "), labels))
      list(unname(unlist(out)))
    }
  }

  if (facet)
    pl <- pl + ggplot2::facet_wrap(c("dx1_name", "dx2_name"),
                labeller = function(labels) label_value2(labels, multi_line = FALSE)) +
                ggplot2::theme(strip.text.x = element_text(size = font_size, angle = 0)) +
                ggplot2::theme(axis.text.x = element_text(size = axis_font_size, angle = x_axis_font_angle)) +
                ggplot2::theme(axis.text.y = element_text(size = axis_font_size))
  # pl <- pl + ggplot2::facet_wrap(~ dx2_name + dx1_name)


  return(pl)
}