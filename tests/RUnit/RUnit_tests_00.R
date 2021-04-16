### --- Test setup ---

if(FALSE) {
  # to automatically enter browser when error shows up:
  # options(error=recover)
  # library("RUnit")
  # library("testthat")
  # library("roxygen2")
  setwd(".."); setwd(".."); getwd()
  library("devtools")

  document()
  load_all("./") # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # stremr:::debug_set() # SET TO DEBUG MODE

  setwd("..");
  install("stremr", build_vignettes = FALSE, dependencies = FALSE) # INSTALL W/ devtools:
  library("stremr")
  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf stremr")  # just create the pdf manual from help files
  # CHECK AND BUILD PACKAGE:
  getwd()
  # setwd("./stremr"); setwd(".."); getwd()
  devtools::check() # runs full check
  devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs faster
  devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  devtools::build(args = "--compact-vignettes") # build package tarball compacting vignettes

  # check reverse dependencies:
  devtools::revdep(dependencies = c("Depends", "Imports", "Suggests", "LinkingTo"),
                    recursive = FALSE, ignore = NULL)
  res <- devtools::revdep_check()
  devtools::revdep_check_summary(res)
  # revdep_check_save_logs(res)

  setwd("..")

  # system("R CMD check --as-cran stremr_0.5.0.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes stremr") # check without building the pdf manual and not building vignettes
      ## system("R CMD build stremr --no-build-vignettes")
      ## system("R CMD build stremr")
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./stremr_0.2.2.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(stremr)
  # stremr:::addvectorfcn("poisson")
  # stremr:::debug_set() # SET TO DEBUG MODE
  # stremr:::debug_off() # SET DEBUG MODE OFF

  # To install a specific branch:
  # options(stremr.verbose = FALSE)
  # devtools::install_github('osofr/gridisl', dependencies = FALSE)
  # devtools::install_github('osofr/simcausal', build_vignettes = FALSE, dependencies = FALSE)
  # devtools::install_github('osofr/stremr', build_vignettes = FALSE, dependencies = FALSE)

  ## devtools version before turbo was renamed to ..turbo in fwrite:
  # devtools::install_github('Rdatatable/data.table', dependencies = FALSE, ref = "4a1f4ba")

}

psi_RDs_DAG2a <- NULL
psi_RDs_DAG2b <- NULL

sample_checks <- function() {   # doesnt run, this is just to show what test functions can be used
  print("Starting tests...")
  checkTrue(1 < 2, "check1")     ## passes fine
  ## checkTrue(1 > 2, "check2")  ## appears as failure in the test protocol
  v <- 1:3
  w <- 1:3
  checkEquals(v, w)               ## passes fine
  names(v) <- c("A", "B", "C")
  ## checkEquals(v, w)            ## fails because v and w have different names
  checkEqualsNumeric(v, w)        ## passes fine because names are ignored
  x <- rep(1:12, 2)
  y <- rep(0:1, 12)
  res <- list(a=1:3, b=letters, LM=lm(y ~ x))
  res2 <- list(a=seq(1,3,by=1), b=letters, LM=lm(y ~ x))
  checkEquals( res, res2)        ## passes fine
  checkIdentical( res, res)
  checkIdentical( res2, res2)
  ## checkIdentical( res, res2)  ## fails because element 'a' differs in type
  fun <- function(x) {
   if(x)
   {
    stop("stop conditions signaled")
   }
   return()
  }
  checkException(fun(TRUE))      ## passes fine
  ## checkException(fun(FALSE))  ## failure, because fun raises no error
  checkException(fun(TRUE), silent=TRUE)
  ##  special constants
  ##  same behaviour as for underlying base functions
  checkEquals(NA, NA)
  checkEquals(NaN, NaN)
  checkEquals(Inf, Inf)
  checkIdentical(NA, NA)
  checkIdentical(NaN, NaN)
  checkIdentical(-Inf, -Inf)
}
`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))

test.GenericModelPredict <- function() {
  options(stremr.verbose = TRUE)

  require("data.table")
  data(OdataNoCENS)
  OdataNoCENS <- as.data.table(OdataNoCENS, key=c("ID", "t"))
  # define lagged N, first value is always 1 (always monitored at the first time point):
  OdataNoCENS[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]
  OdataNoCENS[, "catA" := sample.int(3, size = nrow(OdataNoCENS), replace = TRUE)]

  # -----------------------------------------------------
  # testing fitting and prediction for binary exposure:
  # -----------------------------------------------------
  OData <- importData(OdataNoCENS, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "C", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  gform_CENS <- "C + TI + N ~ highA1c + lastNat1"
  gform_TRT = "TI ~ CVD + highA1c + N.tminus1"
  stratify_TRT = list(TI = c("t == 0L", "t > 0"))
  gform_MONITOR <- "N ~ 1"
  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)
  probA1.class <- OData$modelfit.gA$predict(OData)

  # -----------------------------------------------------
  # same for catogorical exposure:
  # -----------------------------------------------------
  OData <- importData(OdataNoCENS, ID = "ID", t = "t", covars = c("highA1c", "lastNat1"), CENS = "catA", TRT = "TI", MONITOR = "N", OUTCOME = "Y.tplus1")
  gform_TRT = "catA ~ CVD + highA1c + N.tminus1"
  stratify_TRT = list(catA = c("t == 0L", "t > 0"))
  OData <- fitPropensity(OData = OData, gform_CENS = gform_CENS, gform_TRT = gform_TRT, stratify_TRT = stratify_TRT, gform_MONITOR = gform_MONITOR)
  probA1.class <- OData$modelfit.gA$predict(OData)

  options(stremr.verbose = FALSE)
}

test.options <- function() {
  stremrOptions()
  checkException(stremrOptions("blahblah"))
  checkException(stremrOptions("blahblah", 5))
  checkException(stremrOptions("estimator", NULL))
  checkException(stremrOptions("estimator", "blahblah"))
  checkException(stremrOptions("bin_method", "blah"))

  print_stremr_opts()
  stremr:::.onLoad()

  stremr:::gvars$misval
  stremr:::testmisfun()
  stremr:::get.misval()
  stremr:::set.misval(stremr:::gvars, -999)
  stremr:::testmisfun()
  stremr:::get.misval()
  stremr:::set.misval(stremr:::gvars, "-999")
  stremr:::testmisfun()
  stremr:::get.misval()

  stremr:::set.misval(stremr:::gvars, NA)
  stremr:::testmisfun()
  stremr:::get.misval()
}


# test various regression / subsetting schemes and make sure it works as expected
test.regressionCases <- function() {
  # ------------------------------------------------------------------------------------------------------------------
  # 1A) No stratification, TRT=c("A1","A2"), CENS=c("C1","C2","C3"), MONITOR="N":
  # ------------------------------------------------------------------------------------------------------------------
  # Two equivalent ways to specify regression models:
  gform_TRT1 <- c("A1 + A2 ~ L1 + L2")
  gform_TRT2 <- c("A1 ~ L1 + L2", "A2 ~ L1 + L2 + A1")
  # Two equivalent way to specify regressions models:
  gform_CENS1 <- c("C1 + C2 + C3 ~ L1 + L2")
  gform_CENS2 <- c("C1 ~ L1 + L2", "C2 ~ L1 + L2 + C1", "C3 ~ L1 + L2 + C1 + C2")

  gform_MONITOR <- c("N ~ L1 + L2 + A1 + A2")

  # ........
  # stremr <- function(data = Odata, ID = "ID", t = "t",
  #                             covars, CENS = c("C1", "C2", "C3"), TRT = c("A1", "A2"), MONITOR = "N", OUTCOME = "Y",
  #                             gform.CENS = gform_CENS1, gform.TRT = gform_TRT1, gform_MONITOR = gform_MONITOR,
  #                             stratify.CENS = NULL, stratify.TRT = NULL, stratify_MONITOR = NULL)

  # ........
  # stremr <- function(data = Odata, ID = "ID", t = "t",
  #                             covars, CENS = c("C1", "C2", "C3"), TRT = c("A1", "A2"), MONITOR = "N", OUTCOME = "Y",
  #                             gform.CENS = gform_CENS2, gform.TRT = gform_TRT2, gform_MONITOR = gform_MONITOR,
  #                             stratify.CENS = NULL, stratify.TRT = NULL, stratify_MONITOR = NULL)

  # ------------------------------------------------------------------------------------------------------------------
  # 1B) Categorical TRT="A" that codes the same (A1,A2,A3), categorical CENS="C" that codes (C1,C2,C3)
  # Should be equivalent to approach in 1A
  # ------------------------------------------------------------------------------------------------------------------
  gform.TRT3 <- c("A ~ L1 + L2")
  gform.CENS3 <- c("C ~ L1 + L2")
  # ........
  # stremr <- function(data = Odata, ID = "ID", t = "t",
  #                             covars, CENS = c("C"), TRT = c("A"), MONITOR = "N", OUTCOME = "Y",
  #                             gform.CENS = gform.CENS3, gform.TRT = gform.TRT3, gform_MONITOR = gform_MONITOR,
  #                             stratify.CENS = NULL, stratify.TRT = NULL, stratify_MONITOR = NULL)

  # ------------------------------------------------------------------------------------------------------------------
  # 2A) Stratification (subsetting) by logical expressions on TRT=c("A1","A2") & CENS=c("C1","C2","C3"):
  # ******* WHEN ONLY ONE REGRESSION IS SPECIFIED (gform_TRT1 and gform_CENS1) ********
  # ------------------------------------------------------------------------------------------------------------------
  gform_TRT1 <- c("A1 + A2 ~ L1 + L2"); stratify_TRT1 = c("t == 0L", "t > 0L")
  gform_CENS1 <- c("C1 + C2 + C3 ~ L1 + L2"); stratify_CENS1 = c("t == 0L", "t > 0L")
  # ........
  # stremr <- function(data = Odata, ID = "ID", t = "t",
  #                             covars, CENS = c("C1", "C2", "C3"), TRT = c("A1", "A2"), MONITOR = "N", OUTCOME = "Y",
  #                             gform.CENS = gform_CENS1, gform.TRT = gform_TRT1, gform_MONITOR = gform_MONITOR,
  #                             stratify.CENS = stratify_CENS1, stratify.TRT = stratify_TRT1, stratify_MONITOR = NULL)

  # ------------------------------------------------------------------------------------------------------------------
  # 2B) Stratification
  # ******* WHEN > ONE REGRESSION IS SPECIFIED, i.e., gform_TRT2 and gform_CENS2) ********
  # ONE SOLUTION IS TO FORCE stratify.CENS & stratify.TRT TO BE A NAMED LIST.
  # EACH ITEM IN A LIST CORRESPONDS TO VECTOR OF STRATIFICATION RULES
  # ------------------------------------------------------------------------------------------------------------------
  gform_TRT2 <- c("A1 ~ L1 + L2", "A2 ~ L1 + L2 + A1")
  stratify.TRT2 = list(A1=c("t == 0L", "t > 0L"), A2=c("t == 0L", "t > 0L"))

  gform_CENS2 <- c("C1 ~ L1 + L2", "C2 ~ L1 + L2 + C1", "C3 ~ L1 + L2 + C1 + C2")
  stratify_CENS2 = list(C1=c("t == 0L", "t > 0L"), C2=c("t == 0L", "t > 0L"))
  # ........
  # stremr <- function(data = Odata, ID = "ID", t = "t",
  #                             covars, CENS = c("C1", "C2", "C3"), TRT = c("A1", "A2"), MONITOR = "N", OUTCOME = "Y",
  #                             gform.CENS = gform_CENS2, gform.TRT = gform_TRT2, gform_MONITOR = gform_MONITOR,
  #                             stratify.CENS = stratify_CENS2, stratify.TRT = stratify.TRT2, stratify_MONITOR = NULL)

  # ------------------------------------------------------------------------------------------------------------------
  # 2C) Stratification (subsetting) by logical expressions on cat TRT=c("A") and cat CENS="C":
  # ------------------------------------------------------------------------------------------------------------------
  gform.TRT3 <- c("A ~ L1 + L2")
  stratify.TRT3 = c("t == 0L", "t > 0L & A.tminus1 == 0L", "t > 0L & A.tminus1 == 1L")

  gform.CENS3 <- c("C ~ L1 + L2")
  stratify.CENS3 = c("t == 0L", "t > 0L & A.tminus1 == 0L", "t > 0L & A.tminus1 == 1L")
  # ........
  # stremr <- function(data = Odata, ID = "ID", t = "t",
  #                             covars, CENS = "C", TRT = "A", MONITOR = "N", OUTCOME = "Y",
  #                             gform.CENS = gform.CENS3, gform.TRT = gform.TRT3, gform_MONITOR = gform_MONITOR,
  #                             stratify.CENS = stratify.CENS3, stratify.TRT = stratify.TRT3, stratify_MONITOR = NULL)

}


# test that all new nodes added after time-var nodes have a t argument:
test.t.error <- function() {
# ....
}

# testing n.test arg to set.DAG(), including when n.test=0L
test.Nsamp.n.test <- function() {
# ....
}

