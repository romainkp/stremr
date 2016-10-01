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
#' @param TMLE Optional list of a resulting calls to \code{fitTMLE} or a result of a single call to \code{\link{fitTMLE}}.
#' @param GCOMP Optional list of a resulting calls to \code{fitSeqGcomp} or a result of a single call to \code{\link{fitSeqGcomp}}.
#' @param wts_data Optional list of data.tables or a single data.table with weights by regimen.
#' @param SurvByRegimen ... Not implemented ...
#' @param WTtables Table(s) with distribution(s) of the IPTW weights, a result of calling the function \code{\link{get_wtsummary}}
#' @param AddFUPtables Logical, set to \code{TRUE} to print tables describing the distribution of the maximum follow-up times
#' by rule (monitoring and treatment).
#' @param MSM.RDtables List of tables with risk differences returned by the function \code{\link{get_MSM_RDs}}.
#' @param TMLE.RDtables List of tables with risk differences returned by the function \code{\link{get_TMLE_RDs}}.
#' @param format Choose the Pandoc output format for the report file (html, pdf or word).
#' Note that the html report file is always produced in addition to any other selected format.
#' @param skip.modelfits Do not report any of the modeling stats.
#' @param file.name File name for the report file without extension. Default file name is assigned based on the current date.
#' @param file.path Directory path where the report file(s) should be written. Default is to use the system temporary directory.
#' @param openFile Open the report file with OS default viewer?
#' @param keep_md Keep the source .md files?
#' @param keep_tex Keep the source .tex files for pdf output?
#' @param ... Additional arguments may specify the report title (\code{author}), author (\code{title}).
#' Specifying the logical flag \code{only.coefs=TRUE} disables printing of all h2o-specific model summaries.
#' Additional set of arguments control the survival plotting, these are passed on to the function \code{f_plot_survest}:
#' \code{t_int_sel}, \code{y_lab}, \code{x_lab}, \code{miny}, \code{x_legend}, \code{y_legend}.
#' @return String specifying the path to the main report file.
#' @export
make_report_rmd <- function(OData, MSM, NPMSM, TMLE, GCOMP, wts_data, SurvByRegimen,
                            WTtables = NULL, AddFUPtables = FALSE, MSM.RDtables, TMLE.RDtables,
                            format = c("html", "pdf", "word"), skip.modelfits = FALSE,
                            file.name = getOption('stremr.file.name'), file.path = getOption('stremr.file.path'),
                            openFile = TRUE, keep_md = FALSE, keep_tex = FALSE, ...) {
  optArgReport <- list(...)

  if (!rmarkdown::pandoc_available(version = "1.12.3"))
    stop(
"Report functionality requires pandoc (version 1.12.3 or higher).
Please install it.
For more information, go to: http://pandoc.org/installing.html",
call. = FALSE)

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
  # MODEL FITS:
  # -------------------------------------------------------------------------------------
  fitted.coefs.gC <- OData$modelfit.gC$get.fits()
  fitted.coefs.gA <- OData$modelfit.gA$get.fits()
  fitted.coefs.gN <- OData$modelfit.gN$get.fits()

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

  message("writing report to directory: " %+% getwd())
  figure.dir <- file.path(getwd(), "figure/stremr-")
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

  if (openFile) openFileInOS(outfile)

  # resetting directory and other options
  options(opts.bak)
  setwd(wd.bak)
  return(file.path(file.path, outfile))
}
