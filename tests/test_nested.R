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
                                       num.inner.folds=10, max.iters=3,
                                       verbose=FALSE)
  expect_equal(nest.cat[[1]]$panel, c("bmi.cat", "sex.cat"))
  expect_equal(is.na(nest.cat[[1]]$fs$coef), c(TRUE, TRUE, FALSE))
})

context("nested.glm")

test_that("nested.glm",
{
  glm.gauss <- nested.glm(X, Y, c("age", "sex", "bmi", "map"),
                          gaussian(), folds)
  expect_equal(glm.gauss, glm.gauss.ok)

  glm.binom <- nested.glm(X, y.bin, ~ age + sex + bmi + map,
                          binomial(), folds)
  expect_equal(glm.binom, glm.binom.ok)
})
