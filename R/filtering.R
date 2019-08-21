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


#' Compare p(x|y=0) and p(x|y=1) using a Kolmogorov-Smirnov test.
#'
#' @noRd
ks.pval <- function(x, y) {
  is.discrete <- apply(x, 1, function(z) length(unique(z)) < 10)

  ## add some jitter to continuous variables to avoid warnings about ties
  ## note that because of this, results may be slightly different according
  ## to the number of continuous variables in the dataset
  x[!is.discrete, ] <- jitter(x[!is.discrete, ], factor=1e-5)
  x.ctrls <- x[, which(y == 0)]
  x.cases <- x[, which(y == 1)]
  p.values <- matrix(NA, nrow(x), 1)
  for (i in 1:nrow(x)) {
    if (is.discrete[i])
      p.values[i] <- dgof::ks.test(x.ctrls[i, ], ecdf(x.cases[i, ]))$p.value
    else
      p.values[i] <- dgof::ks.test(x.ctrls[i, ], x.cases[i, ])$p.value
  }
  rownames(p.values) <- rownames(x)
  return(-log(p.values))
}

#' Filtering of predictors
#'
#' Filter the predictors, retaining only the top `n`.
#'
#' This performs a univariate test of association of each predictor (not listed
#' in `ignore`) with the outcome, and retains the top `n` variables
#' with smallest p-value according to a Kolmogorov-Smirnov test.
#'
#' @param x Design matrix.
#' @param y Outcome variable.
#' @param n Number of variables to retain.
#' @param ignore Names of variables that should be ignored by the filter.
#'
#' @return
#' An array of indices of variables that are retained by the filter.
#'
#' @keywords internal
filter.predictors <- function(x, y, n, ignore=NULL) {
  stopifnot(n <= ncol(x))
  if (!is.matrix(x))
    x <- sapply(x, as.numeric)
  scores <- ks.pval(t(x[, !colnames(x) %in% ignore]), y)
  limit <- sort(scores, decreasing=TRUE)[n]
  keep.idx <- match(rownames(scores)[scores >= limit], colnames(x))
  return(keep.idx)
}
