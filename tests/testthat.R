Sys.setenv("R_TESTS" = "")
library(testthat)
library(mockery)
test_check("stremr")
