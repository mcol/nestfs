context("arguments")

data(diabetes)
y.gauss <- diabetes$Y
y.binom <- as.integer(diabetes$Y > 140)
y.na <- y.binom
y.na[50] <- NA

test_that("argument checks", {

  expect_error(forward.selection(diabetes, y.binom))
  expect_error(forward.selection(diabetes, y.binom[1:10]))
  expect_error(forward.selection(diabetes, y.na))
  expect_error(forward.selection(diabetes, y.binom, "nonexisting"))

  expect_error(forward.selection(diabetes, y.binom, "age", choose.from=0))
  expect_error(forward.selection(diabetes, y.binom, "age", choose.from=120))

  expect_error(forward.selection(diabetes, y.gauss, family="binomial"))
  expect_error(forward.selection(diabetes, y.gauss, family="poisson"))

  expect_error(forward.selection(diabetes, y.binom, "age", num.inner.folds=1))
  expect_error(forward.selection(diabetes, y.binom, "age", max.iters=0))
  expect_error(forward.selection(diabetes, y.binom, "age", max.pval=0))
  expect_error(forward.selection(diabetes, y.binom, "age", max.pval=1))
  expect_error(forward.selection(diabetes, y.binom, "age", min.llk.diff=-1))
})

context("family validation")
test_that("invalid family inputs",
{
  expect_error(validate.family())
  expect_error(validate.family(NULL))
  expect_error(validate.family(diabetes))
  expect_error(validate.family("nonexisting"))
  expect_error(validate.family(poisson))
})

test_that("valid family inputs",
{
  expect_equal(validate.family("binomial")$family, "binomial")
  expect_equal(validate.family("gaussian")$family, "gaussian")
  expect_equal(validate.family(gaussian())$family, "gaussian")
  expect_equal(validate.family(gaussian)$family,   "gaussian")
})
