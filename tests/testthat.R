library(testthat)
library(nestfs)
library(doParallel)
registerDoParallel()
options(cores=10)

test_file("tests/test_arguments.R")
test_file("tests/test_forward_selection.R")
test_file("tests/test_nested.R")
