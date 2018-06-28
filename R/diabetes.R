#' Diabetes data with interaction terms
#'
#' The dataset consists of 442 patients for which a quantitative measure of
#' disease progression is recorded in the outcome variable \code{Y.diab}. The
#' dataframe of predictors \code{X.diab} includes 10 baseline measurements, in
#' addition to 45 interactions and 9 quadratic terms, for a total of 64
#' variables for each patient.
#'
#' @docType data
#' @name diabetes
#' @aliases X.diab Y.diab
#'
#' @format The dataset consists of the following two variables:
#' \itemize{
#'   \item X.diab: A dataframe of predictors transformed to have zero mean and
#'         unit variance.
#'   \item Y.diab: The outcome variable.
#' }
#'
#' @source
#' B. Efron, T. Hastie, I. Johnstone and R. Tibshirani,
#' "Least angle regression (with discussion)",
#' \emph{Ann. Statist.}, 32 (2), 407-499, 2004.
#' \url{http://www.stanford.edu/~hastie/Papers/LARS/data64.txt}
#'
#' @examples data(diabetes, package="nestfs")
#' @keywords datasets
NULL
