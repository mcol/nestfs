context("arguments")

data(diabetes)
diabetes <- X.diab
y.gauss <- Y.diab
y.binom <- as.integer(Y.diab > 140)
y.short <- y.binom[1:10]

test_that("argument checks",
{
  expect_error(forward.selection(diabetes, y.binom, family=binomial()))
  expect_error(forward.selection(diabetes, y.short, "age", binomial()))
  expect_error(forward.selection(diabetes, y.binom, "nonexisting", binomial()))
  expect_error(forward.selection(within(diabetes, sex[20] <- NA),
                                 y.binom, ~ age + sex, binomial()))

  ## tests for formulas in init.model
  expect_error(forward.selection(diabetes, y.binom, y ~ nonexisting, binomial))
  expect_error(forward.selection(diabetes, y.binom, y ~ NA, binomial))
  expect_error(forward.selection(diabetes, y.binom, "~ age _ lll", binomial))

  ## tests for family
  expect_error(forward.selection(diabetes, y.binom, "age"))
  expect_error(forward.selection(diabetes, y.gauss, "age", binomial()))
  expect_error(forward.selection(diabetes, y.gauss, "age", poisson()))

  ## tests for choose.from
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from=0))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from=120))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from=12.5))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from="nonexisting"))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from=NA))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from=c(TRUE, FALSE)))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 choose.from=diabetes))

  ## tests for num.filter
  expect_error(forward.selection(diabetes, y.gauss, "age", gaussian(),
                                 num.filter=10))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 num.filter=ncol(diabetes)))

  ## tests for other arguments
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 num.inner.folds=1))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 max.iters=0))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 max.pval=0))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 max.pval=1))
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 min.llk.diff=-1))

  ## tests for nested forward selection
  expect_error(nested.forward.selection(diabetes, y.binom, "age"))
  expect_error(nested.forward.selection(diabetes, y.binom, "age", binomial()))
})

context("outcome validation")
test_that("invalid family inputs",
{
  y.binom[50] <- NA
  y.large <- sample(1:2, length(y.binom), replace=TRUE)
  y.categ <- as.factor(y.gauss)
  y.strng <- as.character(y.categ)
  y.inval <- data.frame(matrix(seq(length(y.binom) * 2), nrow=2))
  y.dates <- rep(Sys.Date(), length(y.binom))

  expect_error(forward.selection(diabetes, y.binom, "age", binomial()))
  expect_error(forward.selection(diabetes, y.large, "age", binomial()))
  expect_error(forward.selection(diabetes, y.categ, "age", binomial()))
  expect_error(forward.selection(diabetes, y.strng, "age", gaussian()))
  expect_error(forward.selection(diabetes, y.inval, "age", gaussian()))
  expect_error(forward.selection(diabetes, y.dates, "age", gaussian()))
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
