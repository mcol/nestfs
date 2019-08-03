library(testthat)
library(nestfs)
library(doParallel)
registerDoParallel()
options(cores=2)

test_check("nestfs")
