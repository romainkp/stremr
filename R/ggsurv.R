## -----------------------------------------------------------------------------
## Survival plots with ggplot2
## Functions adapted from GGally R package: https://github.com/ggobi/ggally
## Original R code by Edwin Thoen \email{edwinthoen@@gmail.com},
## modified by Oleg Sofrygin \email{oleg.sofrygin@@gmail.com}
## -----------------------------------------------------------------------------

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("cens", "surv", "up", "low"))
}

#' Survival curves with ggplot2
#'
#' This function produces Kaplan-Meier plots using \code{ggplot2}.
#' As a first argument it needs a \code{survfit} object, created by the
#' \code{survival} package.
#'
#' @export
#' @param estimates A list, one item per regime / intervention. Each list item must be a data.frame containing the
#' survival estimates by time for a single regime / intervention.
#' @param CI should a 95% confidence interval be plotted? Defaults to \code{TRUE}.
#' Uses the standard error estimates provided as a separate column of the input data.
#' @param CI.step When \code{TRUE} the CIs will be plotted as a step function (same as main plot type).
#' When \code{FALSE} the CIs are plotted using \code{ggplot2::geom_ribbon}.
#' @param plot.cens mark the censored observations?
#' @param surv.col colour of the survival estimate. Defaults to black for
#'    one stratum, and to the default \code{ggplot2} colours for multiple
#'    strata. Length of vector with colour names should be either 1 or equal
#'    to the number of strata.
#' @param cens.col colour of the points that mark censored observations.
#' @param lty.est linetype of the survival curve(s). Vector length should be
#'    either 1 or equal to the number of strata.
#' @param lty.ci linetype of the bounds that mark the 95\% CI.
#' @param size.est line width of the survival curve
#' @param size.ci line width of the 95\% CI
#' @param cens.size point size of the censoring points
#' @param cens.shape shape of the points that mark censored observations.
#' @param back.white if \code{TRUE} the background will not be the default
#'    grey of \code{ggplot2}, but will be white with borders around the plot. Defaults to \code{TRUE}
#' @param xlab the label of the x-axis.
#' @param ylab the label of the y-axis.
#' @param main the plot label.
#' @param legend.position Either the coordinates of the legend position inside the plot
#' (e.g., (0.9, 0.2)) or
#' the character word denoting the legend orientation with respect to the plot
#' (e.g., "bottom", "right" or "left").
#' @param surv.name The name of the column containing the survival estimates.
#' @param SE.name The name of the column containing the standard errors (SE) for each time-point estimate of survival.
#' @param order.legend Set to \code{TRUE} to order the legend display by final
#' survival time (highest first).
#' @param ... Additional arguments (not used).
#' @return An object of class \code{ggplot}
#' @author Original R code by Edwin Thoen \email{edwinthoen@@gmail.com}, modified by Oleg Sofrygin \email{oleg.sofrygin@@gmail.com}
ggsurv <- function(
  estimates,
  CI         = TRUE,
  CI.step    = TRUE,
  plot.cens  = TRUE,
  surv.col   = 'gg.def',
  cens.col   = 'gg.def',
  lty.est    = 1,
  lty.ci     = 2,
  size.est   = 0.5,
  size.ci    = size.est,
  cens.size  = 2,
  cens.shape = 3,
  back.white = TRUE,
  xlab       = 'Time',
  ylab       = 'Survival',
  main       = '',
  legend.position = "right",
  surv.name = "St." %+% attr(estimates[[1]], "estimator_short"),
  SE.name = "SE." %+% attr(estimates[[1]], "estimator_short"),
  order.legend = TRUE,
  ...
){

  # legend.position = c(0.1, 0.2),
  # strata <- length(s)

  n.grps <- length(estimates)
  # n <- s$strata
  surv.name <- "St." %+% attr(estimates[[1]], "estimator_short")
  SE.name <- "SE." %+% attr(estimates[[1]], "estimator_short")

  # strataEqualNames <- unlist(strsplit(names(s$strata), '='))
  # ugroups <- strataEqualNames[seq(2, 2*n.grps, by = 2)]

  # getlast <- function(x) {
  #   res <- NULL
  #   maxTime <- max(x$time)
  #   for (mo in names(x$strata)) {
  #     sur <- x[mo]$surv
  #     n <- length(sur)
  #     # grab the last survival value
  #     surValue <- sur[n]
  #     if (isTRUE(all.equal(surValue, 0))) {
  #       # if they die, order by percent complete of max observation.
  #       # tie value of 0 if the last person dies at the last time
  #       surTime <- x[mo]$time[n]
  #       surValue <- (surTime / maxTime) - 1
  #     }
  #     res <- append(res, surValue)
  #   }
  #   return(res)
  # }

  # if (isTRUE(order.legend)) {
  #   group_order <- order(getlast(s), decreasing = TRUE)
  #   lastv <- ugroups[group_order]
  #   if (length(surv.col) == length(n)) {
  #     surv.col <- surv.col[group_order]
  #   }

  #   # if (length(cens.col) == length(n)) {
  #   #   cens.col <- cens.col[group_order]
  #   # }

  # } else {
    # lastv <- ugroups
  # }

  # groups <- factor(ugroups, levels = lastv)
  # gr.name <- strataEqualNames[1]

  gr.name <- "regime"
  gr.df <- vector('list', n.grps)

  # n.ind <- cumsum(c(0, n))

  # for (i in 1:n.grps) {
  #   # indI <- (n.ind[i]+1):n.ind[i+1]

  #   gr.df[[i]] <- data.frame(
  #     time  = c(0, s$time[ indI ]),
  #     surv  = c(1, s$surv[ indI ]),
  #     up    = c(1, s$upper[ indI ]),
  #     low   = c(1, s$lower[ indI ]),
  #     cens  = c(0, s$n.censor[ indI ]),
  #     group = rep(groups[i], n[i] + 1)
  #   )
  # }

  for (i in 1:n.grps) {
    surv_dat <- estimates[[i]]
    gr.df[[i]] <- data.table::data.table(
        time  = surv_dat[["time"]],
        # time  = c(0, s$time[ indI ]),
        surv  = surv_dat[[surv.name]],
        # cens  = c(0, s$n.censor[ indI ]),
        group = surv_dat[["rule.name"]]
        # group = rep(groups[i], n[i] + 1)
      )
      if (SE.name %in% names(surv_dat)) {
        gr.df[[i]][, ("up") := surv + 1.96*surv_dat[[SE.name]]]
        gr.df[[i]][, ("low") := surv - 1.96*surv_dat[[SE.name]]]
      } else {
        CI <- FALSE
      }
  }

  dat      <- data.table::rbindlist(gr.df)

  pl <- ggplot2::ggplot(dat, ggplot2::aes(x = time, y = surv, group = group)) +
    ggplot2::geom_step(ggplot2::aes(col = group, lty = group), size = size.est) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(main)

  pl <- if(surv.col[1] != 'gg.def'){
    scaleValues <- if (length(surv.col) == 1) {
      rep(surv.col, n.grps)
    } else{
      surv.col
    }
    pl + ggplot2::scale_colour_manual(name = gr.name, values = scaleValues)

  } else {
    pl + ggplot2::scale_colour_discrete(name = gr.name)
  }

  lineScaleValues <- if (length(lty.est) == 1) {
    rep(lty.est, n.grps)
  } else {
    lty.est
  }

  pl <- pl + ggplot2::scale_linetype_manual(name = gr.name, values = lineScaleValues)

  if(identical(CI,TRUE)) {
    if(length(surv.col) > 1 || length(lty.est) > 1){
      stop('Either surv.col or lty.est should be of length 1 in order to plot 95% CI with multiple n.grps')
    }

    stepLty <- if ((length(surv.col) > 1 | surv.col == 'gg.def')[1]) {
      lty.ci
    } else {
      surv.col
    }

    if (CI.step) {
      pl <- pl +
        ggplot2::geom_step(ggplot2::aes(y = up, lty = group, col = group), lty = stepLty, size = size.ci) +
        ggplot2::geom_step(ggplot2::aes(y = low,lty = group, col = group), lty = stepLty, size = size.ci)
    } else {
      pl <- pl +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = up, fill = group, linetype = group), alpha = 0.3, size = size.ci, lty = stepLty)
    }

  }

  # if (identical(plot.cens, TRUE) ){
  #   dat.cens <- subset(dat, cens != 0)
  #   dat.cens <- subset(dat.cens, group != "PKD")

  #   if (nrow(dat.cens) == 0) {
  #     stop('There are no censored observations')
  #   }
  #   if (length(cens.col) == 1) {

  #     if (identical(cens.col, "gg.def")) {
  #       # match the colors of the lines
  #       pl <- pl + geom_point(
  #         data = dat.cens,
  #         mapping = aes(y = surv, col = group),
  #         shape = cens.shape,
  #         size    = cens.size,
  #         show.legend = FALSE
  #       )
  #     } else {
  #       # supply the raw color value
  #       pl <- pl + geom_point(
  #         data    = dat.cens,
  #         mapping = aes(y = surv),
  #         shape   = cens.shape,
  #         color   = cens.col,
  #         size    = cens.size
  #       )
  #     }

  #   } else if (length(cens.col) > 0) {
  #     # if(!(identical(cens.col,surv.col) || is.null(cens.col))) {
  #     #   warning ("Color scales for survival curves and censored points don't match.\nOnly one color scale can be used. Defaulting to surv.col")
  #     # }


  #     if (! identical(cens.col, "gg.def")) {
  #       if (length(cens.col) != n.grps) {
  #         warning("Color scales for censored points don't match the number of groups. Defaulting to ggplot2 default color scale")
  #         cens.col <- "gg.def"
  #       }
  #     }

  #     if (identical(cens.col, "gg.def")) {
  #       # match the group color value
  #       pl <- pl + geom_point(
  #         data = dat.cens,
  #         mapping = aes(y = surv, col = group),
  #         shape = cens.shape,
  #         show.legend = FALSE,
  #         size = cens.size
  #       )
  #     } else {

  #       # custom colors and maybe custom shape
  #       uniqueGroupVals = levels(dat.cens$group)
  #       if (length(cens.shape) == 1) {
  #         cens.shape = rep(cens.shape, n.grps)
  #       }

  #       if (length(cens.shape) != n.grps) {
  #         warning("The length of the censored shapes does not match the number of groups (or 1). Defaulting shape = 3 (+)")
  #         cens.shape = rep(3, n.grps)
  #       }
  #       for (i in seq_along(uniqueGroupVals)) {
  #         groupVal = uniqueGroupVals[i]
  #         dtGroup <- subset(dat.cens, group == groupVal)
  #         if (nrow(dtGroup) == 0) {
  #           next
  #         }

  #         pl <- pl + geom_point(
  #           data = dtGroup,
  #           mapping = aes(y=surv),
  #           color = I(cens.col[i]),
  #           shape = cens.shape[i],
  #           show.legend = FALSE,
  #           size = cens.size
  #         )

  #       }
  #     }

  #   }
  # }

  if(identical(back.white, TRUE)) pl <- pl + ggplot2::theme_bw()

  pl <- pl + ggplot2::theme(legend.position = legend.position)

  return(pl)
}


