### --- Test setup ---

if(FALSE) {
  library("RUnit")
  library("roxygen2")
  library("devtools")
  setwd(".."); setwd(".."); getwd()
  document()
  load_all("./") # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # estimtr:::debug_set() # SET TO DEBUG MODE

  setwd("..");
  install("estimtr", build_vignettes = FALSE) # INSTALL W/ devtools:
  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf estimtr")  # just create the pdf manual from help files
  # CHECK AND BUILD PACKAGE:
  getwd()
  # setwd("./estimtr"); setwd(".."); getwd()
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

  # system("R CMD check --as-cran estimtr_0.5.0.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes estimtr") # check without building the pdf manual and not building vignettes
      ## system("R CMD build estimtr --no-build-vignettes")
      ## system("R CMD build estimtr")
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./estimtr_0.2.2.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(estimtr)
  # estimtr:::addvectorfcn("poisson")
  # estimtr:::debug_set() # SET TO DEBUG MODE
  # estimtr:::debug_off() # SET DEBUG MODE OFF

  # To install a specific branch:
  # devtools::install_github('osofr/estimtr', ref = "simnet", build_vignettes = FALSE)
  # options(estimtr.verbose = FALSE)
  # devtools::install_github('osofr/estimtr', build_vignettes = FALSE)
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


# test that all new nodes added after time-var nodes have a t argument:
test.t.error <- function() {
# ....
}

# testing n.test arg to set.DAG(), including when n.test=0L
test.Nsamp.n.test <- function() {
# ....
}

