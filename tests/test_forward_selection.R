context("forward selection")
options(cores=10)

data(diabetes)
X <- diabetes[, 2:10]
Y <- diabetes$Y
y.binom <- as.integer(Y > 140)

source("dump_forward_selection.R")
fs.gauss.ok$call <- fs.binom.ok$call <- NULL

test_that("linear regression",
{
  fs.gauss <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                                num.inner.folds=10, max.iters=3)
  fs.gauss$call <- NULL
  expect_equal(fs.gauss, fs.gauss.ok)

  fs.gauss <- forward.selection(X, Y, c("age"), family="gaussian",
                                choose.from=c("sex"))
  expect_is(fs.gauss, "fs")
  expect_equal(fs.gauss$final.model, "age")
  expect_named(fs.gauss$iter1, c("median.diff.llk", "total.diff.llk", "p.value"))
  expect_equal(dim(fs.gauss$all.iter[[1]]), c(2, 30))
  expect_length(fs.gauss$all.iter, 1)
  expect_length(fs.gauss$panel, 0)
})

test_that("logistic regression",
{
  fs.binom <- forward.selection(X, y.binom, c("age", "sex"), family="binomial",
                                num.inner.folds=10, max.iters=3)
  fs.binom$call <- NULL
  expect_equal(fs.binom, fs.binom.ok)

  ## logical outcome variable
  y.binom <- diabetes$Y > 140
  fs.logic <- forward.selection(X, y.binom, c("age", "sex"), family="binomial",
                                num.inner.folds=10, max.iters=3)
  fs.logic$call <- NULL
  expect_equal(fs.logic, fs.binom.ok)

  ## factor outcome variable
  y.binom <- factor(diabetes$Y > 140)
  fs.logic <- forward.selection(X, y.binom, c("age", "sex"), family="binomial",
                                num.inner.folds=10, max.iters=3)
  fs.logic$call <- NULL
  expect_equal(fs.logic, fs.binom.ok)
})

test_that("init.model",
{
  ## empty init.model
  fs.1 <- forward.selection(X, Y, c(), family="gaussian",
                            num.inner.folds=10, max.iters=3)
  fs.2 <- forward.selection(X, Y, y ~ 1, family="gaussian",
                            num.inner.folds=10, max.iters=3)
  fs.1$call <- fs.2$call <- NULL
  expect_equal(fs.1, fs.2)
  expect_equal(fs.1$panel, c("bmi", "ltg", "map"))

  ## formula as a character string
  fs.1 <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                            num.inner.folds=10, max.iters=1)
  fs.2 <- forward.selection(X, Y, "y ~ age + sex", family="gaussian",
                            num.inner.folds=10, max.iters=1)
  fs.1$call <- fs.2$call <- NULL
  expect_equal(fs.1, fs.2)

  ## formula with misspecified or no left-hand side
  fs.3 <- forward.selection(X, Y, Z ~ age + sex, family="gaussian",
                            num.inner.folds=10, max.iters=1)
  fs.4 <- forward.selection(X, Y, ~ age + sex, family="gaussian",
                            num.inner.folds=10, max.iters=1)
  fs.5 <- forward.selection(X, Y, "~ age + sex", family="gaussian",
                            num.inner.folds=10, max.iters=1)
  fs.3$call <- fs.4$call <- fs.5$call <- NULL
  expect_equal(fs.1, fs.3)
  expect_equal(fs.1, fs.4)
  expect_equal(fs.1, fs.5)

  ## models with interaction terms
  fs.1 <- forward.selection(X, Y, y ~ age*sex, family="gaussian",
                            num.inner.folds=10, max.iters=3)
  fs.2 <- forward.selection(X, Y, "age*sex", family="gaussian",
                            num.inner.folds=10, max.iters=3)
  fs.3 <- forward.selection(X, Y, c("age", "sex", "age:sex"), family="gaussian",
                            num.inner.folds=10, max.iters=3)
  fs.1$call <- fs.2$call <- fs.3$call <- NULL
  expect_equal(fs.1, fs.2)
  expect_equal(fs.1, fs.3)
})

test_that("choose.from",
{
  ## no new variables to choose from
  fs.gauss <- forward.selection(X, Y, c("age"), family="gaussian",
                                choose.from=c("age"))
  expect_length(fs.gauss$panel, 0)
  expect_length(fs.gauss$iter1, 0)
  expect_length(fs.gauss$all.iter, 0)

  ## choose.from with indices and variable names
  fs.1 <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                            choose.from=c(3:6),
                            num.inner.folds=10, max.iters=3)
  fs.2 <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                            choose.from=c("bmi", "map", "tc", "ldl"),
                            num.inner.folds=10, max.iters=3)
  fs.1$call <- fs.2$call <- NULL
  expect_equal(fs.1, fs.2)
})

context("nested forward selection")

folds <- create.folds(3, nrow(X))

test_that("nested forward selection",
{
  nest.1 <- nested.forward.selection(X, Y, ~ age + sex, gaussian(), folds,
                                     num.inner.folds=10, max.iters=3)
  nest.2 <- nested.forward.selection(X, Y, ~ age + sex, gaussian(), folds[1],
                                     num.inner.folds=10, max.iters=3)
  nest.1[[1]]$call <- NULL
  nest.2[[1]]$call <- NULL
  expect_equal(nest.1[[1]], nest.2[[1]])
})
