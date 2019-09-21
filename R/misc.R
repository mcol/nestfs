##=============================================================================
##
## Copyright (c) 2018-2019 Marco Colombo
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##=============================================================================


#' Validate model formula and dataset and extract the outcome
#'
#' @param model Formula to be checked.
#' @param data Data frame or matrix containing the variables used in the model.
#'
#' @return
#' A numeric vector containing the outcome variable.
#'
#' @noRd
validate.model.outcome <- function(model, data) {
  if (is.character(model) && length(model) > 1)
    stop("Model formula specified incorrectly.", call.=FALSE)
  model <- as.formula(model)
  tt <- terms(model)
  if (attr(tt, "response") == 0)
    stop("No outcome variable specified in the model.", call.=FALSE)
  if (!inherits(data, c("data.frame", "matrix")))
    stop("'data' must be a data frame or a matrix.", call.=FALSE)
  if (any(!all.vars(model) %in% colnames(data)))
    stop("Not all model variables are in 'data'.", call.=FALSE)
  mf <- model.frame(model, data)
  y <- validate.outcome(model.response(mf))
  return(y)
}

#' Validate the outcome variable
#'
#' Ensure that the outcome variable has been specified correctly.
#'
#' @param y Outcome variable to test.
#'
#' @return
#' A valid outcome variable. The function throws an error if the outcome
#' variable cannot be used.
#'
#' @noRd
validate.outcome <- function(y) {
  if (anyNA(y))
    stop("Outcome variable contains missing values.", call.=FALSE)
  if (is.character(y))
    stop("Outcome variable cannot be a character vector.", call.=FALSE)
  if (is.factor(y)) {
    if (nlevels(y) != 2)
      stop("A factor outcome variable can only have two levels.", call.=FALSE)
    y <- as.integer(y) - 1
  }
  if (!(is.numeric(y) || is.logical(y)))
    stop("Outcome variable of invalid type.", call.=FALSE)

  return(as.numeric(y))
}

#' Validate initial model
#'
#' Ensure that the initial model has been specified correctly.
#'
#' @param model Model definition to test.
#'
#' @return
#' A formula describing the initial model. The function throws an error if the
#' model parameter cannot be used.
#'
#' @importFrom methods is
#' @noRd
validate.init.model <- function(model) {
  if (is.null(model) || length(model) == 0) {
    model <- y ~ 1
  }
  else if (is.character(model)) {
    if (any(model == ""))
      stop("init.model contains an empty string.", call.=FALSE)
    if (length(model) == 1 && grepl("~", model))
      model <- as.formula(model)
    else
      model <- reformulate(model, "y")
  }
  else if (!is(model, "formula"))
    stop("init.model specified incorrectly.", call.=FALSE)

  ## rename the left-hand side or add it if not present
  model <- update(model, "nestfs_y_ ~ .")

  return(model)
}

#' Validate the family argument
#'
#' Ensure that the family argument has been specified correctly.
#' This is inspired by code in `glm`.
#'
#' @param family Family argument to test.
#' @param y Outcome variable.
#'
#' @return
#' A valid family. The function throws an error if the family argument cannot
#' be used.
#'
#' @importFrom methods is
#' @noRd
validate.family <- function(family, y) {
  if (missing(family))
    stop("Argument of 'family' is missing.", call.=FALSE)
  if (is.character(family))
    tryCatch(
      family <- get(family, mode="function", envir=parent.frame(2)),
      error=function(e)
        stop("'", family, "' is not a valid family.", call.=FALSE)
    )
  if (is.function(family))
    family <- family()
  if (!is(family, "family"))
    stop("Argument of 'family' is not a valid family.", call.=FALSE)
  if (!family$family %in% c("gaussian", "binomial"))
    stop("Only 'gaussian' and 'binomial' are supported families.", call.=FALSE)

  if (family$family == "gaussian" && is.factor(y))
    stop("Factor outcome variable not valid with family=gaussian().", call.=FALSE)
  if (family$family == "binomial") {
    if (length(table(y)) != 2)
      stop("y must contain two classes with family=binomial().", call.=FALSE)
    if (!is.factor(y) && any(y < 0 | y > 1))
      stop("y must contain 0-1 values with family=binomial().", call.=FALSE)
  }

  return(family)
}

#' Validate the choose.from argument
#'
#' Ensure that the `choose.from` argument has been specified correctly.
#'
#' @param choose.from Argument to test.
#' @param x Dataframe of predictors.
#'
#' @return
#' A valid vector of variable indices. The function throws an error if the
#' argument cannot be used.
#'
#' @noRd
validate.choose.from <- function(choose.from, x) {
  if (is.null(choose.from))
    choose.from <- seq(ncol(x))
  else {
    if (is.numeric(choose.from)) {
      if (anyNA(choose.from))
        stop("choose.from contains missing values.", call.=FALSE)
      if (length(choose.from) > 0 &&
          (min(choose.from) < 1 || max(choose.from) > ncol(x)))
        stop("choose.from contains out of bounds indices.", call.=FALSE)
      if (any(choose.from != as.integer(choose.from)))
        stop("choose.from contains floating point values.", call.=FALSE)
    }
    else if (is.character(choose.from)) {
      choose.from <- match(choose.from, colnames(x))
      if (anyNA(choose.from))
        stop("choose.from contains names that cannot be matched.", call.=FALSE)
    }
    else
      stop("choose.from should be an integer or character vector.", call.=FALSE)
  }
  return(as.integer(choose.from))
}

#' Validate the folds argument
#'
#' Ensure that the `folds` argument has been specified correctly.
#'
#' @param folds Argument to test.
#'
#' @return
#' A valid list of folds. The function throws an error if the argument cannot
#' be used.
#'
#' @noRd
validate.folds <- function(folds, x) {
  if (!is.list(folds))
    stop("folds expected to be a list.", call.=FALSE)
  all.idx <- unlist(folds)
  if (anyNA(all.idx))
    stop("folds contains missing values.", call.=FALSE)
  if (!is.numeric(all.idx))
    stop("folds contains non-numerical values.", call.=FALSE)
  if (any(all.idx != as.integer(all.idx)))
    stop("folds contains non-integer values", call.=FALSE)
  if (any(table(all.idx) > 1))
    stop("folds contains repeated indices.", call.=FALSE)
  if (any(all.idx <= 0 | all.idx > nrow(x)))
    stop("folds contains out of bounds indices.", call.=FALSE)
  return(folds)
}
