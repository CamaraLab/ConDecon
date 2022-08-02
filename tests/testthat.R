Sys.setenv("R_TESTS" = "")
library(testthat)
library(ConDecon)

test_check("ConDecon")
