data(diabetes)
X <- X.diab[, 1:9]
Y <- Y.diab
y.bin <- as.integer(Y > 140)
folds <- create.folds(2, nrow(X), 0)

source("dump_nested.R")

context("nested forward selection")

test_that("nested forward selection",
{
  nest.gauss <- nested.forward.selection(X, Y, ~ age + sex, gaussian, folds,
                                         num.inner.folds=10, max.iters=3,
                                         verbose=FALSE)
  expect_equal(nest.gauss, nest.gauss.ok)

  nest.binom <- nested.forward.selection(X, y.bin, ~ age + sex, binomial, folds,
                                         num.inner.folds=10, max.iters=3,
                                         verbose=FALSE)
  expect_equal(nest.binom, nest.binom.ok)
})

context("nested.glm")

test_that("nested.glm",
{
  glm.gauss <- nested.glm(X[, c("age", "sex", "bmi", "map")], Y,
                          gaussian(), folds)
  expect_equal(glm.gauss, glm.gauss.ok)

  glm.binom <- nested.glm(X[, c("age", "sex", "bmi", "map")], y.bin,
                          binomial(), folds)
  expect_equal(glm.binom, glm.binom.ok)
})
