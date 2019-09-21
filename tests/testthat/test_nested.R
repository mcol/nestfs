data(diabetes)
X <- diabetes[, 1:9]
Y <- diabetes$Y
XY <- cbind(X, Y)
y.bin <- factor(Y > 140)
folds <- create.folds(2, nrow(X), 0)

source("dump_nested.R")

context("nested forward selection")

test_that("nested forward selection",
{
  nest.gauss <- nested.fs(Y ~ age + sex, XY, gaussian, folds,
                          num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(nest.gauss, nest.gauss.ok)

  summ.gauss <- summary(nest.gauss)
  expect_named(summ.gauss, c("vars", "percent", "coef", "coefIQR",
                             "rank", "rankIQR", "diffLogLik", "diffLogLikIQR"))
  expect_equivalent(summ.gauss$rank, c(2, 1, 2, 3, 3))
  expect_equivalent(summ.gauss$coefIQR[1], "(28.74, 30.17)")
  expect_equivalent(summ.gauss$diffLogLikIQR[1], "(24.38, 40.36)")

  summ.gauss.iter1 <- summary(nest.gauss, iter1=TRUE)
  expect_named(summ.gauss.iter1,
               c("med.diff.llk", "iqr.diff.llk", "med.pvalue", "iqr.pvalue"))
  expect_equal(rownames(summ.gauss.iter1),
               c("bmi", "ltg", "tch", "map", "hdl", "tc", "ldl"))
  expect_equivalent(summ.gauss.iter1$med.pvalue,
                    c(0.001, 0.001, 0.004, 0.003, 0.013, 0.213, 0.199))

  nest.binom <- nested.fs(y.bin ~ age + sex, cbind(X, y.bin), binomial, folds,
                          num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(nest.binom, nest.binom.ok)

  ## interaction terms present in the initial model
  nest.int <- nested.fs(Y ~ age * sex, XY, gaussian, folds,
                        num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_true("age:sex" %in% rownames(nest.int[[1]]$model))

  ## test the selection of categorical variables
  X$sex.cat <- factor(X$sex)
  X$bmi.cat <- factor(cut(X$bmi, c(-Inf, -1, 1, Inf)))
  nest.cat <- nested.fs(Y ~ age, cbind(X, Y), gaussian, folds,
                        choose.from=c("sex.cat", "bmi.cat"), min.llk.diff=0,
                        num.inner.folds=10, max.iters=3, verbose=FALSE)
  expect_equal(nest.cat[[1]]$panel, c("bmi.cat"))
  expect_equal(is.na(nest.cat[[1]]$fs$coef), c(TRUE, TRUE))

  summ.cat <- summary(nest.cat)
  expect_equivalent(summ.cat$coef, c(NA, -15.711))
  expect_equivalent(summ.cat$diffLogLik, c(25.599, 0.009))
})

test_that("nested forward selection old interface",
{
  nest.gauss <- nested.forward.selection(X, Y, ~ age + sex, gaussian, folds,
                                         num.inner.folds=10, max.iters=3,
                                         verbose=FALSE)
  expect_equal(nest.gauss, nest.gauss.ok)

  summ.gauss <- summary(nest.gauss)
  expect_named(summ.gauss, c("vars", "percent", "coef", "coefIQR",
                             "rank", "rankIQR", "diffLogLik", "diffLogLikIQR"))
  expect_equivalent(summ.gauss$rank, c(2, 1, 2, 3, 3))
  expect_equivalent(summ.gauss$coefIQR[1], "(28.74, 30.17)")
  expect_equivalent(summ.gauss$diffLogLikIQR[1], "(24.38, 40.36)")

  summ.gauss.iter1 <- summary(nest.gauss, iter1=TRUE)
  expect_named(summ.gauss.iter1,
               c("med.diff.llk", "iqr.diff.llk", "med.pvalue", "iqr.pvalue"))
  expect_equal(rownames(summ.gauss.iter1),
               c("bmi", "ltg", "tch", "map", "hdl", "tc", "ldl"))
  expect_equivalent(summ.gauss.iter1$med.pvalue,
                    c(0.001, 0.001, 0.004, 0.003, 0.013, 0.213, 0.199))

  nest.binom <- nested.forward.selection(X, y.bin, ~ age + sex, binomial, folds,
                                         num.inner.folds=10, max.iters=3,
                                         verbose=FALSE)
  expect_equal(nest.binom, nest.binom.ok)

  ## interaction terms present in the initial model
  nest.int <- nested.forward.selection(X, Y, ~ age * sex, gaussian, folds,
                                       num.inner.folds=10, max.iters=3,
                                       verbose=FALSE)
  expect_true("age:sex" %in% rownames(nest.int[[1]]$model))

  ## test the selection of categorical variables
  X$sex.cat <- factor(X$sex)
  X$bmi.cat <- factor(cut(X$bmi, c(-Inf, -1, 1, Inf)))
  nest.cat <- nested.forward.selection(X, Y, ~ age, gaussian, folds,
                                       choose.from=c("sex.cat", "bmi.cat"),
                                       num.inner.folds=10, max.iters=3, min.llk.diff=0,
                                       verbose=FALSE)
  expect_equal(nest.cat[[1]]$panel, c("bmi.cat"))
  expect_equal(is.na(nest.cat[[1]]$fs$coef), c(TRUE, TRUE))

  summ.cat <- summary(nest.cat)
  expect_equivalent(summ.cat$coef, c(NA, -15.711))
  expect_equivalent(summ.cat$diffLogLik, c(25.599, 0.009))
})

context("nested.glm")

test_that("nested.glm",
{
  glm.gauss <- nested.glm(Y ~ age + sex + bmi + map,
                          diabetes, gaussian(), folds)
  expect_equal(glm.gauss, glm.gauss.ok)

  glm.binom <- nested.glm(y.bin ~ age + sex + bmi + map, cbind(diabetes, y.bin),
                          binomial(), folds)
  expect_equal(glm.binom, glm.binom.ok)

  ## performance
  perf.gauss <- nested.performance(glm.gauss)
  expect_named(perf.gauss, c("observed", "predicted", "performance"))
  expect_named(attributes(perf.gauss), c("names", "measure", "class"))
  expect_s3_class(perf.gauss, "nestperf")
  expect_equal(perf.gauss$performance, 0.62174173)
  expect_length(perf.gauss$observed, length(Y))
  expect_output(print(perf.gauss), "Correlation coefficient:")

  perf.binom <- nested.performance(glm.binom)
  expect_equal(perf.binom$performance, 0.78288733)
  expect_output(print(perf.binom), "Area under the curve:")
})
