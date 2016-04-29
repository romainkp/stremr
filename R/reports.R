#' @import ggplot2
#' @import plyr
# @import pander
# @importFrom pander evalsOptions
NULL

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


#' @export
make_report_rmd <- function(OData, MSM, MSM.list, Surv.byregimen, format = "html", file.name = getOption('stremr.file.name'), file.path = getOption('stremr.file.path')) {
  sVartypes <- gvars$sVartypes

  if (!missing(MSM)) {
    # if ()
    # handle separately if MSM is a list of many MSMs -> will need to do the plotting for each MSM in the list
    Surv.byregimen <- MSM$St
  }

  if (!missing(MSM.list)) {
    # need to figure out how to handle it -> define some variables which will control report-scrit-rmd.R
    # ...
  }

  # -------------------------------------------------------------------------------------
  # DISCRIPTIVE STATISTICS:
  # -------------------------------------------------------------------------------------
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
  # NEED TO ADD *****NUMBER OF OBSERVATIONS USED**** for each model within each:
  fitted.coefs.gC <- OData$modelfit.gC$get.fits(format_table = TRUE)
  fitted.coefs.gA <- OData$modelfit.gA$get.fits(format_table = TRUE)
  fitted.coefs.gN <- OData$modelfit.gN$get.fits(format_table = TRUE)

  # outvar = self$outvar, predvars = self$predvars, stratify = self$subset_expr)
  # -------------------------------------------------------------------------------------
  # **** NEED TO ADD RD tables ****
  # -------------------------------------------------------------------------------------
  RD.IPAW_tperiod1 <- MSM.IPAW$RD.IPAW_tperiod1
  RD.IPAW_tperiod2 <- MSM.IPAW$RD.IPAW_tperiod2
  RR.IPAW_tperiod1 <- MSM.IPAW$RR.IPAW_tperiod1
  RR.IPAW_tperiod2 <- MSM.IPAW$RR.IPAW_tperiod2

  ## path issue on Windows
  file.path     <- gsub('\\', '/', file.path, fixed = TRUE)
  # find the full path to the report template:
  report.file <- system.file('report', "report-script-rmd.R", package = 'stremr')

  ## set working directory where to write the report:
  opts.bak <- options()                      # backup options
  wd.bak   <- getwd()
  setwd(file.path)

  format_pandoc <- format %+% "_document"
  outfile <- file.name %+% "." %+% ifelse(format %in% "word", "docx", format)

  print("writing report to directory: " %+% getwd())

  figure.dir <- file.path(getwd(), "figure/stremr-")
  # output_format = "html_document",
  report.html <- tryCatch(rmarkdown::render(report.file,
                          output_dir = getwd(), intermediates_dir = getwd(), output_file = file.name%+%".html", clean = TRUE,
                          output_options = list(keep_md = TRUE, toc = TRUE, toc_float = TRUE,
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
        output_options <- list(keep_tex = TRUE, toc = TRUE, number_sections = TRUE, fig_caption = TRUE,
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

  openFileInOS(outfile)
  # resetting directory and other options
  options(opts.bak)
  setwd(wd.bak)
}



#' @export
make_report <- function(OData, Surv.byregimen, file.name = getOption('stremr.file.name'), file.path = getOption('stremr.file.path')) {
  sVartypes <- gvars$sVartypes
  ## path issue on Windows
  file.path     <- gsub('\\', '/', file.path, fixed = TRUE)

  # find the full path to the report template:
  report.file <- system.file('report', "graphs.brew", package = 'stremr')
  # read in the report file:
  # txt <- readLines(report.file, warn = FALSE, encoding = 'UTF-8') # load template from file path


  ## pregenerate file name
  # if (grepl('%T', file.name)) file.name <- gsub('%T', gsub('\\\\|/|:|\\.', '-', fp), file.name, fixed = TRUE)
  # file.name <- gsub('--', '-', file.name, fixed = TRUE)
  # if (grepl('%N', file.name)) {
  #   if (length(strsplit(sprintf('placeholder%splaceholder', file.name), '%N')[[1]]) > 2) stop('File name contains more then 1 "%N"!')
  #   similar.files <-  list.files(file.path(file.path, 'plots'), pattern = sprintf('^%s\\.(jpeg|tiff|png|svg|bmp)$', gsub('%t', '[a-z0-9]*', gsub('%N|%n|%i', '[[:digit:]]*', file.name))))
  #   if (length(similar.files) > 0) {
  #     similar.files <- sub('\\.(jpeg|tiff|png|svg|bmp)$', '', similar.files)
  #     rep <- gsub('%t|%n|%i', '[a-z0-9]*', strsplit(basename(file.name), '%N')[[1]])
  #     `%N` <- max(as.numeric(gsub(paste(rep, collapse = '|'), '', similar.files))) + 1
  #   } else
  #     `%N` <- 1
  #   file.name <- gsub('%N', `%N`, file.name, fixed = TRUE)
  # }

  ## set working directory where to write the report:
  opts.bak <- options()                      # backup options
  wd.bak   <- getwd()
  setwd(file.path)

  # browser()
  # evalsOptions('graph.name', file.name)
  # assign('.rapport.body', paste(b, collapse = '\n'), envir = e)
  # assign('.graph.name', file.name, envir = e)
  # assign('.graph.dir', evalsOptions('graph.dir'), envir = e)
  # assign('.graph.hi.res', graph.hi.res, envir = e)
  # if (grepl("w|W", .Platform$OS.type)) # we are on Windows
  #   assign('.tmpout', 'NUL', envir = e)
  # else
  #   assign('.tmpout', '/dev/null', envir = e)
  # report <- tryCatch(eval(parse(text = 'Pandoc.brew(text = .rapport.body, graph.name = .graph.name, graph.dir = .graph.dir, graph.hi.res = .graph.hi.res, output = .tmpout)'), envir = e), error = function(e) e)

  ## Initialize a new Pandoc object
  myReport <- Pandoc$new()

  ## Add author, title and date of document
  myReport$author <- 'Gergely Daróczi'
  myReport$title  <- 'Demo'
  ## Add some free text
  myReport$add.paragraph('Hello there, this is a really short tutorial!')
  ## Add maybe a header for later stuff
  myReport$add.paragraph('# Showing some raw R objects below')

  f_plot_survest <- function(surv_res_est, t_int_sel, y_lab, miny) {
    ptsize <- 0.7
    counter <- 0
    if (missing(y_lab)) y_lab <- ""
    if (missing(t_int_sel)) t_int_sel <- seq_along(surv_res_est[[1]])
    if (missing(miny)) miny <- min(unlist(lapply(surv_res_est, min)))

    for(d.j in names(surv_res_est)){
      counter <- counter+1
      plot(t_int_sel,surv_res_est[[d.j]][t_int_sel],col=counter,type='b',cex=ptsize,ylim=c(miny,1),
        ylab = y_lab, xlab="Quarter since study entry")
      par(new=TRUE)
    }
    legend(12,0.96,legend=names(surv_res_est),col=c(1:length(names(surv_res_est))), cex=ptsize, pch=1)
  }
  # graph.dir = 'my_plots',
  # evals('f_plot_survest(Surv.byregimen)', graph.output = 'png')[[1]]$result
  # myReport$add.paragraph(evals('f_plot_survest(Surv.byregimen)', graph.output = 'png')[[1]]$result)
  evalsOptions('graph.unify', FALSE)
  # browser()
  # str(myReport)
  myReport$add(f_plot_survest(Surv.byregimen))
  print(file.path)
  print(file.name)
  # print(myReport)
  myReport$format <- 'md'
  myReport$export(file.name, open = FALSE)

  # report.file <- system.file('report', "graphs.brew", package = 'stremr')
  report.file <- system.file('report', "short-code-long-report.brew", package = 'stremr')

  Pandoc.brew(file = report.file, output = file.name %+% "." %+% myReport$format, convert = FALSE,
              open = FALSE,
              envir = parent.frame(), append = TRUE)
              # graph.name, graph.dir, graph.hi.res = FALSE, text = NULL,

  Pandoc.convert(file.name %+% "." %+% myReport$format, format = 'html', open = TRUE)
  getwd()
  # resetting directory and other options
  options(opts.bak)
  setwd(wd.bak)
}
