context("forward selection")
options(cores=10)

data(diabetes)
X <- diabetes[, 2:10]
Y <- diabetes$Y
y.binom <- as.integer(Y > 140)

source("dump_forward_selection.R")

test_that("linear regression",
{
  fs.gauss <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                                num.inner.folds=10, max.iters=3)
  expect_equal(fs.gauss, fs.gauss.ok)
})

test_that("logistic regression",
{
  fs.binom <- forward.selection(X, y.binom, c("age", "sex"), family="binomial",
                                num.inner.folds=10, max.iters=3)
  expect_equal(fs.binom, fs.binom.ok)
})
