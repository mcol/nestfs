##=============================================================================
##
## Copyright (c) 2014-2019 Marco Colombo
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
