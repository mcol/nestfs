context("forward selection")

data(diabetes)
X <- X.diab[, 1:9]
Y <- Y.diab
y.binom <- as.integer(Y > 140)

source("dump_forward_selection.R")

test_that("linear regression",
{
  fs.gauss <- forward.selection(X, Y, c("age", "sex"), family="gaussian",
                                num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.gauss, fs.gauss.ok)

  fs.gauss <- forward.selection(X, Y, c("age"), family="gaussian",
                                choose.from=c("sex"), verbose=FALSE)
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
                                num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.binom, fs.binom.ok)
})

test_that("init.model",
{
  ## empty init.model
  fs.1 <- forward.selection(X, Y, c(), family="gaussian",
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1$init, character(0))
  expect_equal(fs.1$panel, c("bmi", "ltg", "map"))
  expect_equal(fs.1$init.model, "1")
  expect_equal(fs.1$final.model, "bmi + ltg + map")

  ## models with interaction terms
  fs.1 <- forward.selection(X, Y, y ~ age*sex, family="gaussian",
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1$init, c("age", "sex"))
  expect_equal(fs.1$panel, c("bmi", "ltg", "map"))
  expect_equal(fs.1$init.model, "age + sex + age:sex")
  expect_equal(fs.1$final.model, "age + sex + bmi + ltg + map + age:sex")

  ## summary
  summ <- summary(fs.1)
  expect_named(summ, c("vars", "fdr", "llks", "diffs", "iter"))
  expect_equal(summ$llks[1:4], c(NA, -2136.283195, -2047.874234, -2004.104712))
  expect_equal(summ$iter, c(NA, NA, 1:3))
})

test_that("choose.from",
{
  ## no new variables to choose from
  fs.gauss <- forward.selection(X, Y, c("age"), family="gaussian",
                                choose.from=c("age"), verbose=FALSE)
  expect_length(fs.gauss$panel, 0)
  expect_length(fs.gauss$iter1, 0)
  expect_length(fs.gauss$all.iter, 0)

  ## choose.from with indices and variable names
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
  fs.1 <- forward.selection(X, y.binom, ~ age + sex, binomial(),
                            num.filter=5,
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.1$fs, fs.binom.ok$fs)
  expect_equal(dim(fs.1$iter1), c(5, 3))

  fs.2 <- forward.selection(X.diab, y.binom, ~ age + sex, binomial(),
                            choose.from=30:40, num.filter=5,
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_length(fs.2$panel, 0)
  expect_equal(dim(fs.2$iter1), c(5, 3))

  fs.3 <- forward.selection(X.diab, y.binom, ~ age + sex, binomial(),
                            num.filter=5, filter.ignore=c("bmi", "nonexisting"),
                            num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(fs.3$fs, fs.binom.ok$fs)
  expect_equal(dim(fs.3$iter1), c(6, 3))
})
