### --- Test setup ---
if(FALSE) {
  # CHECK AND BUILD PACKAGE:
  # library("stremr")
  library("roxygen2")
  library("devtools")
  library("testthat")
  library("data.table")
  setwd(".."); setwd(".."); getwd()
  document()
  load_all("./", create = FALSE) # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # stremr:::debug_set() # SET TO DEBUG MODE

  setwd("..");
  install("stremr", build_vignettes = FALSE, dependencies = FALSE) # INSTALL W/ devtools:

  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf stremr")  # just create the pdf manual from help files

  getwd()
  # setwd("./stremr"); setwd(".."); getwd()
  devtools::check(args = "--as-cran")
  devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs faster
  # devtools::check() # runs check with devtools
  # devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs check with devtools
  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  # devtools::build()
  devtools::build(args = "--compact-vignettes")
  # devtools::build(args = c("--compact-vignettes", "--resave-data"))
  devtools::build_win(args = "--as-cran") # build package on CRAN servers (windows os?)
  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  # devtools::build(args = "--compact-vignettes") # build package tarball compacting vignettes
  # devtools::build(args = "--no-build-vignettes") # build package tarball compacting vignettes
  # devtools::build() # build package tarball
  setwd("..")
  # system("R CMD check --as-cran stremr_0.1.0.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes stremr") # check without building the pdf manual and not building vignettes
      ## system("R CMD build stremr --no-build-vignettes --as-cran")
      # system("R CMD build stremr --resave-data")
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./stremr_0.2.0.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(stremr)
  # stremr:::debug_set() # SET TO DEBUG MODE
  # stremr:::debug_off() # SET DEBUG MODE OFF

  # To install a specific branch:
  # devtools::install_github('osofr/simcausal', ref = "simnet", build_vignettes = FALSE)
  # devtools::install_github('osofr/stremr', ref = "master", build_vignettes = FALSE)

  # TEST COVERATE:
  # if your working directory is in the packages base directory
  # package_coverage()
  # or a package in another directory
  # cov <- package_coverage("stremr")
  # view results as a data.frame
  # as.data.frame(cov)
  # zero_coverage() can be used to filter only uncovered lines.
  # zero_coverage(cov)
}


