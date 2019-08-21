##=============================================================================
##
## Copyright (c) 2013-2019 Marco Colombo
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


#' Cross-validation folds
#'
#' Create a list of indices corresponding to cross-validation folds.
#'
#' @param num.folds Number of folds to be created.
#' @param num.rows Number of observations in the dataset.
#' @param seed Seed of the random number generator. If `NULL`, the folds
#'        generated will be different at each invocation; for reproducibility
#'        of results, it is recommended to set this to a specific value.
#'
#' @return
#' A list of length `num.folds` containing the indices of the observations
#' to be withdrawn for validation in each fold.
#'
#' @note
#' Note that the number of observations withdrawn in each fold may not be
#' exactly the same if `num.folds` is not an integer divisor of `num.rows`.
#'
#' @examples
#' all.folds <- create.folds(50, 307, 0)
#' @export
create.folds <- function(num.folds, num.rows, seed=NULL) {
  if (!is.null(seed))
    set.seed(seed)
  ## ensure that the test indices are not repeated
  full.idx <- numeric(ceiling(num.rows / num.folds) * num.folds)
  full.idx[1:num.rows] <- sample(num.rows)
  fold.matrix <- as.data.frame(matrix(full.idx, ncol=num.folds, byrow=TRUE))
  folds <- lapply(fold.matrix, function(z) sort(z[z > 0]))
  names(folds) <- NULL
  return(folds)
}
