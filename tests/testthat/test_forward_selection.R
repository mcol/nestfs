context("forward selection")

data(diabetes)
X <- diabetes[, 1:9]
Y <- diabetes$Y
XY <- cbind(X, Y)
y.binom <- as.integer(Y > 140)

source("dump_forward_selection.R")

test_that("linear regression",
{
  fs.gauss <- fs(Y ~ age + sex, XY, family="gaussian",
                 num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.gauss, fs.gauss.ok)

  fs.gauss <- fs(Y ~ age, diabetes, family="gaussian",
                 choose.from="sex", verbose=FALSE)
  expect_is(fs.gauss, "fs")
  expect_equal(fs.gauss$final.model, "age")
  expect_named(fs.gauss$iter1, c("median.diff.llk", "total.diff.llk", "p.value"))
  expect_equal(dim(fs.gauss$all.iter[[1]]), c(2, 30))
  expect_length(fs.gauss$all.iter, 1)
  expect_length(fs.gauss$panel, 0)
})

test_that("logistic regression",
{
  fs.binom <- fs(y.binom ~ age + sex, cbind(X, y.binom), family="binomial",
                 num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.binom, fs.binom.ok)
})

test_that("init.model",
{
  ## empty init.model
  fs.1 <- fs(Y ~ 1, XY, family="gaussian",
             num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- forward.selection(X, Y, c(), family="gaussian",
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)
  expect_equal(fs.1$init, character(0))
  expect_equal(fs.1$panel, c("ltg", "bmi", "map"))
  expect_equal(fs.1$init.model, "1")
  expect_equal(fs.1$final.model, "ltg + bmi + map")

  ## models with interaction terms
  fs.1 <- fs(Y ~ age*sex, XY, family="gaussian",
             num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- forward.selection(X, Y, y ~ age*sex, family="gaussian",
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)
  expect_equal(fs.1$init, c("age", "sex"))
  expect_equal(fs.1$panel, c("ltg", "bmi", "map"))
  expect_equal(fs.1$init.model, "age + sex + age:sex")
  expect_equal(fs.1$final.model, "age + sex + ltg + bmi + map + age:sex")

  ## summary
  summ <- summary(fs.1)
  expect_named(summ, c("vars", "fdr", "llks", "diffs", "iter"))
  expect_equal(summ$llks[1:4], c(NA, -2137.599696, -2056.442947, -2004.431740))
  expect_equal(summ$iter, c(NA, NA, 1:3))
})

test_that("choose.from",
{
  ## no new variables to choose from
  fs.gauss <- fs(Y ~ age, diabetes, family="gaussian",
                 choose.from="age", verbose=FALSE)
  fs.gauss2 <- forward.selection(X, Y, c("age"), family="gaussian",
                                 choose.from=c("age"), verbose=FALSE)
  expect_equal(fs.gauss, fs.gauss2)
  expect_length(fs.gauss$panel, 0)
  expect_length(fs.gauss$iter1, 0)
  expect_length(fs.gauss$all.iter, 0)

  ## choose.from with indices and variable names
  fs.1 <- fs(Y ~ age + sex, diabetes, family="gaussian",
             choose.from=5, num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- fs(Y ~ age + sex, diabetes, family="gaussian",
             choose.from="tc", num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)

  fs.1 <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                            choose.from=5,
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                            choose.from="tc",
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)
})

test_that("univariate filter",
{
  fs.1 <- fs(y.binom ~ age + sex, cbind(X, y.binom), binomial(),
             num.filter=5, num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- forward.selection(X, y.binom, ~ age + sex, binomial(),
                            num.filter=5,
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)
  expect_equal(fs.1$fs, fs.binom.ok$fs)
  expect_equal(dim(fs.1$iter1), c(5, 3))

  fs.1 <- fs(y.binom ~ age + sex, cbind(diabetes, y.binom), binomial(),
             choose.from=30:40, num.filter=5, min.llk.diff=0,
             num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- forward.selection(diabetes, y.binom, ~ age + sex, binomial(),
                            choose.from=30:40, num.filter=5, min.llk.diff=0,
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)
  expect_length(fs.1$panel, 3)
  expect_equal(dim(fs.1$iter1), c(5, 3))

  diab <- diabetes[, -match("Y", colnames(diabetes))]
  fs.1 <- fs(y.binom ~ age + sex, cbind(diab, y.binom), binomial(),
             num.filter=5, filter.ignore=c("bmi", "nonexisting"),
             num.inner.folds=10, max.iters=3, verbose=FALSE)
  fs.2 <- forward.selection(diab, y.binom, ~ age + sex, binomial(),
                            num.filter=5, filter.ignore=c("bmi", "nonexisting"),
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1, fs.2)
  expect_equal(fs.1$fs, fs.binom.ok$fs)
  expect_equal(dim(fs.1$iter1), c(6, 3))
})
