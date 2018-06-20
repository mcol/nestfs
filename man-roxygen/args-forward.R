#' @param x.all Dataframe of predictors: this should include all variables in
#'        the initial set and the variables that are allowed to enter the
#'        selected panel.
#' @param y.all Outcome variable: if \code{family="binomial"}, it is expected
#'        to only have \code{0-1} entries.
#' @param init.vars Initial set of variables (ignored if \code{init.model} is
#'        not \code{NULL}).
