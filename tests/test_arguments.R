context("arguments")

data(diabetes)
diabetes <- X.diab
y.gauss <- Y.diab
y.binom <- as.integer(Y.diab > 140)
y.short <- y.binom[1:10]

test_that("argument checks",
{
  expect_error(forward.selection(diabetes, y.binom, family=binomial()),
               "is missing, with no default")
  expect_error(forward.selection(diabetes, y.short, "age", binomial()),
               "Mismatched dimensions")
  expect_error(forward.selection(diabetes, y.binom, "nonexisting", binomial()),
               "not present in x")
  expect_error(forward.selection(within(diabetes, sex[20] <- NA),
                                 y.binom, ~ age + sex, binomial()),
               "Missing values in the variables")

  ## tests for init.model
  expect_error(forward.selection(diabetes, y.binom, y ~ nonexisting, binomial),
               "not present in x")

  ## tests for family
  expect_error(forward.selection(diabetes, y.binom, "age"),
               "Argument of 'family' is missing")

  ## tests for num.filter
  expect_error(forward.selection(diabetes, y.gauss, "age", gaussian(),
                                 num.filter=10),
               "only be used with family=binomial()")
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 num.filter=ncol(diabetes)),
               "cannot exceed the number of available predictors")

  ## tests for filter.ignore
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 num.filter=5, filter.ignore=NA),
               "should be a character vector or NULL")
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 num.filter=5, filter.ignore=1:10),
               "should be a character vector or NULL")

  ## tests for other arguments
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 num.inner.folds=1),
               "should be at least 5")
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 max.iters=0),
               "should be at least 1")
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 max.pval=0),
               "should be between 0 and 1")
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 max.pval=1),
               "should be between 0 and 1")
  expect_error(forward.selection(diabetes, y.binom, "age", binomial(),
                                 min.llk.diff=-1),
               "cannot be negative")

  ## tests for nested forward selection
  expect_error(nested.forward.selection(diabetes, y.binom, "age"),
               "is missing, with no default")
  expect_error(nested.forward.selection(diabetes, y.binom, "age", binomial()),
               "is missing, with no default")
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

  expect_error(forward.selection(diabetes, y.binom, "age", binomial()),
               "contains missing values")
  expect_error(forward.selection(diabetes, y.categ, "age", binomial()),
               "can only have two levels")
  expect_error(forward.selection(diabetes, y.strng, "age", gaussian()),
               "cannot be a character vector")
  expect_error(forward.selection(diabetes, y.inval, "age", gaussian()),
               "invalid type")
  expect_error(forward.selection(diabetes, y.dates, "age", gaussian()),
               "invalid type")
})

context("init.model validation")
test_that("invalid init.model inputs",
{
  expect_error(validate.init.model(NA),
               "specified incorrectly")
  expect_error(validate.init.model(y ~ NA),
               "invalid model formula in ExtractVars")
  expect_error(validate.init.model("~ NA"),
               "invalid model formula in ExtractVars")
  expect_error(validate.init.model("~ age _ lll"),
               "unexpected input")
  expect_error(validate.init.model(c("age", "*sex")),
               "unexpected '*'")
  expect_error(validate.init.model(c("age", "")),
               "contains an empty string")
})

test_that("valid init.model inputs",
{
  expect_equal(validate.init.model(NULL),
               nestfs_y_ ~ 1)
  expect_equal(validate.init.model(character(0)),
               nestfs_y_ ~ 1)
  expect_equal(validate.init.model(integer(0)),
               nestfs_y_ ~ 1)
  expect_equal(validate.init.model(~ age + sex),
               nestfs_y_ ~ age + sex)
  expect_equal(validate.init.model(y ~ age * sex),
               nestfs_y_ ~ age + sex + age:sex)
  expect_equal(validate.init.model(c("age", "sex")),
               nestfs_y_ ~ age + sex)
  expect_equal(validate.init.model("age*sex"),
               nestfs_y_ ~ age + sex + age:sex)
  expect_equal(validate.init.model("age:sex"),
               nestfs_y_ ~ age:sex)
  expect_equal(validate.init.model("z ~ age * sex"),
               nestfs_y_ ~ age + sex + age:sex)
  expect_equal(validate.init.model(~ notexisting),
               nestfs_y_ ~ notexisting)
})

context("family validation")
test_that("invalid family inputs",
{
  expect_error(validate.family(),
               "Argument of 'family' is missing")
  expect_error(validate.family(NULL),
               "is not a valid family")
  expect_error(validate.family(diabetes),
               "is not a valid family")
  expect_error(validate.family("nonexisting"),
               "is not a valid family")
  expect_error(validate.family(poisson),
               "are supported families")
  expect_error(validate.family(binomial(), y.gauss),
               "must contain two classes")
  expect_error(validate.family(binomial(), y.binom + 1),
               "must contain 0-1 values")
  expect_error(validate.family(binomial(), y.binom - 1),
               "must contain 0-1 values")
})

test_that("valid family inputs",
{
  expect_equal(validate.family("binomial", c(1, 0))$family, "binomial")
  expect_equal(validate.family("gaussian")$family, "gaussian")
  expect_equal(validate.family(gaussian())$family, "gaussian")
  expect_equal(validate.family(gaussian)$family,   "gaussian")
})

context("choose.from validation")
test_that("invalid choose.from inputs",
{
  expect_error(validate.choose.from(0, diabetes),
               "out of bounds indices")
  expect_error(validate.choose.from(120, diabetes),
               "out of bounds indices")
  expect_error(validate.choose.from(12.5, diabetes),
               "contains floating point values")
  expect_error(validate.choose.from("nonexisting", diabetes),
               "names that cannot be matched")
  expect_error(validate.choose.from(c(1:5, NA), diabetes),
               "contains missing values")
  expect_error(validate.choose.from(NA, diabetes),
               "integer or character vector")
  expect_error(validate.choose.from(c(TRUE, FALSE), diabetes),
               "integer or character vector")
  expect_error(validate.choose.from(diabetes, diabetes),
               "integer or character vector")
})

test_that("valid choose.from inputs",
{
  expect_equal(validate.choose.from(NULL, diabetes),
               1:ncol(diabetes))
  expect_equal(validate.choose.from(character(0), diabetes),
               integer(0))
  expect_equal(validate.choose.from(integer(0), diabetes),
               integer(0))
  expect_equal(validate.choose.from(colnames(diabetes)[c(1, 6, 7)], diabetes),
               c(1, 6, 7))
  expect_equal(validate.choose.from(5,  diabetes),
               validate.choose.from(5L, diabetes))
})

context("folds validation")
test_that("invalid folds inputs",
{
  expect_error(validate.folds(NULL),
               "expected to be a list")
  expect_error(validate.folds(1:100),
               "expected to be a list")
  expect_error(validate.folds(list(c(1:5), NA)),
               "contains missing values")
  expect_error(validate.folds(list(c(1:5), c("a", "b"))),
               "contains non-numerical values")
  expect_error(validate.folds(list(c(1:5), c(3.2, 4.5))),
               "contains non-integer values")
  expect_error(validate.folds(list(c(1:5), c(5:10))),
               "contains repeated indices")
  expect_error(validate.folds(list(c(1:5), c(6:10, 0)), diabetes),
               "contains out of bounds indices")
  expect_error(validate.folds(list(c(1:5), c(6, nrow(diabetes) + 1)), diabetes),
               "contains out of bounds indices")
})
