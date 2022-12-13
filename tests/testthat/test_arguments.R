data(diabetes)
y.gauss <- diabetes$Y
y.binom <- as.integer(diabetes$Y > 140)
y.short <- y.binom[1:10]
folds <- create.folds(2, nrow(diabetes), 0)

context("arguments")
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
  expect_error(fs(Y ~ nonexisting, diabetes, gaussian()),
               "Not all model variables are in 'data'")
  expect_error(forward.selection(diabetes, y.binom, y ~ nonexisting, binomial),
               "not present in x")
  expect_error(fs(y ~ ., diabetes, gaussian()),
               "No selection possible")
  expect_error(forward.selection(diabetes, y.binom, y.binom ~ ., binomial),
               "No selection possible")
  expect_error(fs("y ~ .", diabetes, gaussian()),
               "No selection possible")
  expect_error(forward.selection(diabetes, y.binom, "y.binom ~ .", binomial),
               "No selection possible")
  expect_error(fs(y ~ age + ., diabetes, gaussian()),
               "No selection possible")
  expect_error(fs(y ~ age + sex, diabetes[, c("age", "sex")], gaussian()),
               "No selection possible")

  ## tests for family
  expect_error(fs(Y ~ age, diabetes),
               "Argument of 'family' is missing")
  expect_error(forward.selection(diabetes, y.binom, "age"),
               "Argument of 'family' is missing")

  ## tests for num.filter
  expect_error(fs(Y ~ age, diabetes, gaussian(), num.filter=10),
               "only be used with family=binomial()")
  expect_error(fs(I(Y > 140) ~ age, diabetes, binomial(),
                  num.filter=ncol(diabetes)),
               "cannot exceed the number of available predictors")
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
               "Argument of 'family' is missing")
  expect_error(nested.forward.selection(diabetes, y.binom, "age", binomial()),
               "is missing, with no default")

  ## tests for nested.glm
  expect_error(nested.glm(Y ~ age),
               "is missing, with no default")
  expect_error(nested.glm(Y ~ age, diabetes),
               "Argument of 'family' is missing")
  expect_error(nested.glm(Y ~ nonexisting, diabetes, gaussian(), folds),
               "Not all model variables are in 'data'")

  ## tests for nested.performance
  expect_error(nested.performance(),
               "is missing, with no default")
  expect_error(nested.performance(NULL),
               "Object is not of 'nestfs' or 'nestglm' class")
})

context("outcome validation")
test_that("invalid outcome inputs",
{
  expect_error(validate.outcome(c(1:5, NA, 6:10)),
               "contains missing values")
  expect_error(validate.outcome(as.factor(y.gauss)),
               "can only have two levels")
  expect_error(validate.outcome(as.character(y.gauss)),
               "cannot be a character vector")
  expect_error(validate.outcome(data.frame(matrix(1:40, nrow=2))),
               "invalid type")
  expect_error(validate.outcome(rep(Sys.Date(), 10)),
               "invalid type")
  expect_error(validate.outcome(NULL),
               "invalid type")
})

test_that("valid outcome inputs",
{
  expect_is(validate.outcome(c(1.1, 2.2, 3.3, 4.4)), "numeric")
  expect_equal(validate.outcome(c(1.1, 2.2, 3.3, 4.4)),
               c(1.1, 2.2, 3.3, 4.4))
  expect_equal(validate.outcome(c(1L, 0L, 0L, 1L)),
               c(1, 0, 0, 1))
  expect_equal(validate.outcome(factor(c("Yes", "No", "No", "Yes"))),
               c(1, 0, 0, 1))
  expect_equal(validate.outcome(c(TRUE, FALSE, FALSE, TRUE)),
               c(1, 0, 0, 1))
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
  expect_error(validate.init.model("y ~ ."),
               "No selection possible")
  expect_error(validate.init.model(c("age", "*sex")),
               "unexpected '*'")
  expect_error(validate.init.model(c("age", "")),
               "contains an empty string")
})

test_that("valid init.model inputs",
{
  expect_is(validate.init.model(NULL), "formula")
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
  expect_error(validate.family(binomial()),
               "is missing, with no default")
  expect_error(validate.family(binomial(), y.gauss),
               "must contain two classes")
  expect_error(validate.family(binomial(), y.binom + 1),
               "must contain 0-1 values")
  expect_error(validate.family(binomial(), y.binom - 1),
               "must contain 0-1 values")
  expect_error(validate.family(gaussian(), as.factor(y.binom)),
               "Factor outcome variable not valid with family=gaussian()")
})

test_that("valid family inputs",
{
  expect_is(validate.family(gaussian(), y.gauss), "family")
  expect_equal(validate.family("binomial", c(1, 0))$family, "binomial")
  expect_equal(validate.family("binomial", factor(c(1, 0)))$family, "binomial")
  expect_equal(validate.family("gaussian", y.gauss)$family, "gaussian")
  expect_equal(validate.family(gaussian(), y.gauss)$family, "gaussian")
  expect_equal(validate.family(gaussian, y.gauss)$family,   "gaussian")
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
  expect_is(validate.choose.from(5, diabetes), "integer")
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

test_that("valid folds inputs",
{
  expect_is(validate.folds(list(1:221, 222:nrow(diabetes)), diabetes), "list")
})
