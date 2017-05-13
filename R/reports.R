# @import pander
# @importFrom pander evalsOptions
NULL

# nocov start
#' Open file
#'
#' Tries to open a file with operating system's default program.
#' @param f file (with full path)
#' @references This function is a fork of David Hajage's \code{convert} function: \url{https://github.com/eusebe/ascii/blob/master/R/export.r}
#' @export
openFileInOS <- function(f) {
  if (missing(f)) {
    stop('No file to open!')
  }
  f <- path.expand(f)
  if (!file.exists(f)) {
    stop('File not found!')
  }
  if (grepl('w|W', .Platform$OS.type)) {
    # we are on Windows
    shell.exec(f) #nolint
  } else {
    if (grepl('darwin', version$os)) {
      # Mac
      system(paste(shQuote('open'), shQuote(f)), wait = FALSE, ignore.stderr = TRUE)
    } else {
      # Linux-like
      system(paste(shQuote('/usr/bin/xdg-open'), shQuote(f)), #nolint
             wait = FALSE,
             ignore.stdout = TRUE)
    }
  }
}
# nocov end

# ---------------------------------------------------------------------------------------------
#' Generate report(s) with modeling stats and survival estimates using pandoc.
#'
#' @param OData Input data object returned by the function \code{\link{importData}}.
#' @param MSM The MSM object fits returned by the function \code{\link{survMSM}}.
#' @param NPMSM Optional list of a resulting calls to \code{survNPMSM} or a result of a single call to \code{\link{survNPMSM}}.
#' @param TMLE Optional list of a resulting calls to \code{fit_TMLE} or a result of a single call to \code{\link{fit_TMLE}}.
#' @param GCOMP Optional list of a resulting calls to \code{fit_GCOMP} or a result of a single call to \code{\link{fit_GCOMP}}.
#' @param WTtables Table(s) with distribution(s) of the IPTW weights, a result of calling the function \code{\link{get_wtsummary}}
#' @param FUPtables Subject-specific \code{data.table} with maximum follow-up time saved for each subject in the column named \code{max.t}.
#' See \code{\link{get_FUPtimes}} for additional details.
#' @param MSM.RDtables List of tables with risk differences returned by the function \code{\link{get_MSM_RDs}}.
#' @param TMLE.RDtables List of tables with risk differences returned by the function \code{\link{get_TMLE_RDs}}.
#' @param plotKM Logical, set to \code{TRUE} to plot KM survival curves when \code{NPMSM} argument is specified. Default is \code{FALSE}.
#' @param printEstimateTables ...
#' @param format Choose the Pandoc output format for the report file (html, pdf or word).
#' Note that the html report file is always produced in addition to any other selected format.
#' @param skip.modelfits Do not report any of the modeling stats.
#' @param file.name File name for the report file without extension. Default file name is assigned based on the current date.
#' @param file.path Directory path where the report file(s) should be written. Default is to use the system temporary directory.
#' @param openFile Open the report file with OS default viewer?
#' @param keep_md Keep the source .md files?
#' @param keep_tex Keep the source .tex files for pdf output?
#' @param serve_html_rmote Serve the html report as a webpage via R package "rmote".
# ' Requires prior initialization of the back-end server with rmote::start_rmote()
#' @param save_report_data Save the data needed for generating this report as a list in a separate 'report_name.Rd' file.
#' @param use_ggplot Set to \code{TRUE} to enable plotting with ggplot.
#' @param ... Additional arguments may specify the report title (\code{author}), author (\code{title}).
#' Specifying the logical flag \code{only.coefs=TRUE} disables printing of all h2o-specific model summaries.
#' Additional set of arguments control the survival plotting, these are passed on to the function \code{f_plot_survest}:
#' \code{t_int_sel}, \code{y_lab}, \code{x_lab}, \code{miny}, \code{x_legend}, \code{y_legend}.
#' @return String specifying the path to the main report file.
#' @export
make_report_rmd <- function(OData, MSM, NPMSM, TMLE, GCOMP,
                            WTtables, FUPtables, MSM.RDtables, TMLE.RDtables,
                            plotKM = FALSE, printEstimateTables = FALSE,
                            format = c("html", "pdf", "word"), skip.modelfits = FALSE,
                            file.name = getOption('stremr.file.name'), file.path = getOption('stremr.file.path'),
                            openFile = TRUE, serve_html_rmote = FALSE, keep_md = FALSE, keep_tex = FALSE, save_report_data = FALSE, use_ggplot = TRUE, ...) {
  optArgReport <- list(...)

  if (!rmarkdown::pandoc_available(version = "1.12.3"))
    stop(
"Report functionality requires pandoc (version 1.12.3 or higher).
Please install it.
For more information, go to: http://pandoc.org/installing.html",
call. = FALSE)

  clean_est_object <- function(est_obj) {
    if ("wts_data" %in% names(est_obj)) {
      est_obj$wts_data <- NULL
    } else if ("wts_data" %in% names(est_obj[[1]])) {
      for (idx in seq_along(est_obj)) est_obj[[idx]]$wts_data <- NULL
    }

    if ("IC.Var.S.d" %in% names(est_obj)) {
      est_obj$IC.Var.S.d <- NULL
    } else if ("IC.Var.S.d" %in% names(est_obj[[1]])) {
      for (idx in seq_along(est_obj)) est_obj[[idx]]$IC.Var.S.d <- NULL
    }
    return(est_obj)
  }

  if ("author" %in% names(optArgReport)) {
    author <- optArgReport[['author']]
    assert_that(is.character(author))
  } else {
    author <- "Insert Author"
  }

  if ("title" %in% names(optArgReport)) {
    title <- optArgReport[['title']]
    assert_that(is.character(author))
  } else {
    title <- "stremr Analysis Report"
  }

  if ("only.coefs" %in% names(optArgReport)) {
    only.coefs <- optArgReport[['only.coefs']]
    assert_that(is.logical(only.coefs))
  } else {
    only.coefs <- FALSE
  }

  ## -------------------------------------------------------------------------------------
  ## MODEL FITS:
  ## -------------------------------------------------------------------------------------
  if (!skip.modelfits) {
    model_fits_gC <- OData$modelfit.gC$get.fits()
    model_summaries_gC <- OData$modelfit.gC$get.model.summaries()

    model_fits_gA <- OData$modelfit.gA$get.fits()
    model_summaries_gA <- OData$modelfit.gA$get.model.summaries()

    model_fits_gN <- OData$modelfit.gN$get.fits()
    model_summaries_gN <- OData$modelfit.gN$get.model.summaries()
  }

  # model_fits_gC[[1]]$get_overall_best_model()[[1]]
  # model_fits_gC[[1]]$get_best_models(K=1)
  # # # try(pander::pander(model_fits_gN[[1]]$getMSEtab, caption = "Overall Performance by Model"))
  # h2o_gridobj <- model_fits_gC[[1]]$get_modelfits_grid()[[1]]
  # capture.output(print(h2o_gridobj))
  # xgb_gridobj <- model_fits_gC[[1]]$get_modelfits_grid()[[3]]
  # # capture.output(print(xgb_gridobj))
  # # # class(xgb_gridobj)
  # # # is.data.frame(xgb_gridobj)
  # xgb_gridobj <- xgb_gridobj[ , names(xgb_gridobj)[!(names(xgb_gridobj) %in% c("glob_params", "xgb_fit", "fit", "params"))], with = FALSE]
  # try(pander::pander_return(xgb_gridobj, caption = "Grid Details"))

  # # class(model_fits_gC[[1]]$get_modelfits_grid()[[1]])
  # # str(model_fits_gC[[1]]$get_modelfits_grid()[[1]])
  # # showMethods(class = "H2OGrid", printTo = FALSE )
  # pander::pander(capture.output(print(h2o_gridobj))[-2])

  # model_fits_gC[[1]]$get_modelfits_grid()
  # model_fits_gC[[1]]$get_best_models()

  # model_fits_gC[[1]]$OData_train
  # model_fits_gC[[1]]$OData_valid
  # grids <- model_fits_gN[[1]]$get_modelfits_grid()

  # best_model <- model_fits_gC[[1]]$get_best_models(K=1)[[1]]
  # str(best_model@model)
  # str(best_model@model$training_metrics)
  # best_model@model$training_metrics
  # gridisl:::pander.H2OBinomialMetrics(best_model@model$training_metrics, "train")
  # gridisl::print_tables(best_model)

  # str(fitted.coefs.gC[[1]])
  # fitted.coefs.gC[[1]]$show(print_format = FALSE, model_stats = TRUE, all_fits = TRUE)
  # fitted.coefs.gC[[1]]$summary(all_fits = TRUE)
  # str(fitted.coefs.gC[[1]]$get_best_models()[[1]])
  # fitted.coefs.gC[[1]]$get_best_models()[[1]]
  # print(fitted.coefs.gC[[1]]$get_best_models(), only.coefs = TRUE)

  ## -------------------------------------------------------------------------------------
  ## Number of unique ID and number of person time obs
  ## -------------------------------------------------------------------------------------
  nuniqueIDs <- OData$nuniqueIDs
  nobs <- OData$nobs
  nuniquets <- OData$nuniquets
  t_name <- OData$nodes$tnode
  use_ggplot <- use_ggplot

  ## -------------------------------------------------------------------------------------
  ## Create report data object (list) to be saved along with the report itself
  ## -------------------------------------------------------------------------------------
  if (save_report_data) {
    OData_save <- OData$clone()
    OData_save$emptydat.sVar
    report_results_list <- list(OData = OData_save)

    if (!missing(MSM)) {
      wts_data <- MSM$wts_data
      MSM$wts_data <- NULL
      MSM$IC.Var.S.d <- NULL
      report_results_list <- c(report_results_list, list(MSM = MSM))
    }

    if (!missing(NPMSM)) {
      NPMSM <- clean_est_object(NPMSM)
      report_results_list <- c(report_results_list, list(NPMSM = NPMSM))
    }

    if (!missing(TMLE)) {
      TMLE <- clean_est_object(TMLE)
      report_results_list <- c(report_results_list, list(TMLE = TMLE))
    }

    if (!missing(GCOMP)) {
      GCOMP <- clean_est_object(GCOMP)
      report_results_list <- c(report_results_list, list(GCOMP = GCOMP))
    }

    if (!missing(WTtables)) report_results_list <- c(report_results_list, list(WTtables = WTtables))
    if (!missing(FUPtables)) report_results_list <- c(report_results_list, list(FUPtables = FUPtables))
    if (!missing(MSM.RDtables)) report_results_list <- c(report_results_list, list(MSM.RDtables = MSM.RDtables))
    if (!missing(TMLE.RDtables)) report_results_list <- c(report_results_list, list(TMLE.RDtables = TMLE.RDtables))
  }

  # -------------------------------------------------------------------------------------
  # TO DO: DISCRIPTIVE STATISTICS
  # -------------------------------------------------------------------------------------
  # *** Hisogram of follow-up times by rule (based on the last value of t for each ID);
  # *** A table of summary(fuptimes), with quantiles, min, max, ...
  # sVartypes <- gvars$sVartypes
  # For each covariate in OData$nodes, create either:
    # 1) a density plot for each  OData$type.sVar[varnames] %in% sVartypes$cont
    # 2) a histogram plot for each variable OData$type.sVar[varnames] %in% c(sVartypes$cont, sVartypes$cat, sVartypes$bin)
    # 3) a frequency table for each c(sVartypes$cat, sVartypes$bin)
  # For each factor coverted to binary dummies from OData$new.factor.names do the same
  # * Report number lost to follow-up by time for nodes$Cnodes, by nodes$tnode
  # * Report number lost to follow-up by nodes$Ynode, by nodes$tnode
  # * Report number in each exposure cat in nodes$Anodes, by nodes$tnode

  # -------------------------------------------------------------------------------------
  # RD tables
  # RDs.IPAW.tperiods <- MSM$RDs.IPAW.tperiods
  # -------------------------------------------------------------------------------------
  ## path issue on Windows
  file.path     <- gsub('\\', '/', file.path, fixed = TRUE)
  # find the full path to the report template:
  report.file <- system.file('report', "report-script-rmd.R", package = 'stremr')

  ## set working directory where to write the report:
  opts.bak <- options() # backup options
  wd.bak   <- getwd()
  setwd(file.path)

  format <- format[1L]
  format_pandoc <- format %+% "_document"
  outfile <- file.name %+% "." %+% ifelse(format %in% "word", "docx", format)
  figure.subdir <- "figure." %+% file.name
  # figure.dir <- file.path(getwd(), "figure/stremr-")

  message("writing report to directory: " %+% getwd())
  message("writing related figures to report subdirectory: " %+% figure.subdir)

  figure.dir <- file.path(getwd(), figure.subdir %+% "/stremr-")
  report.html <- tryCatch(rmarkdown::render(report.file,
                          output_dir = getwd(), intermediates_dir = getwd(), output_file = file.name%+%".html", clean = TRUE,
                          output_options = list(keep_md = keep_md, toc = TRUE, toc_float = TRUE,
                                                number_sections = TRUE, fig_caption = TRUE,
                                                # mathjax = "local", self_contained = FALSE)
                                                md_extensions = "+escaped_line_breaks")
                          ),
                  error = function(e) e)

  if (inherits(report.html, 'error')) {
    options(opts.bak)
    setwd(wd.bak)
    stop(report.html$message)
  }

  if (!format %in% "html"){
      if (format %in% "pdf") {
        output_options <- list(keep_tex = keep_tex, toc = TRUE, number_sections = TRUE, fig_caption = TRUE,
                               md_extensions = "+escaped_line_breaks")
        # output_options <- list(keep_tex = TRUE, pandoc = list(arg="+escaped_line_breaks"))
        # output_options <- list(keep_tex = TRUE, pandoc = list(arg="markdown+escaped_line_breaks"))
      } else {
        output_options <- NULL
      }
    report.other <- tryCatch(rmarkdown::render(report.file,
                            output_dir = getwd(), intermediates_dir = getwd(), output_file = outfile,
                            output_format = format_pandoc, clean = TRUE,
                            output_options = output_options),
                    error = function(e) e)
    if (inherits(report.other, 'error')) {
      options(opts.bak)
      setwd(wd.bak)
      stop(report.other$message)
    }
  }

  if (save_report_data) {
    # print(object.size(report_results_list), units = "MB")
    save(list = "report_results_list", file = file.path(file.path, file.name) %+% ".Rd")
  }

  if (openFile) openFileInOS(outfile)
  if (serve_html_rmote) {
    stop("must install development version of stremr from github for this to work")
  #     reqrmote <- requireNamespace("rmote", quietly = TRUE)
  #     serv_tmle_exists <- exists("serve_rmd_html", where = "package:rmote")
  #     if (!reqrmote || !serv_tmle_exists) stop("Please install the latest version of 'rmote' package by typing this into the terminal:
  # devtools::install_github('hafen/rmote')
  # ", call. = FALSE)
  #     rmote::serve_rmd_html(file.path, file.name%+%".html")
  }

  # resetting directory and other options
  options(opts.bak)
  setwd(wd.bak)
  return(file.path(file.path, outfile))
}
