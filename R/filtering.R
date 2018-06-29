#' Compare p(x|y=0) and p(x|y=1) using a Kolmogorov-Smirnov test.
#'
#' @importFrom dgof ks.test
#' @importFrom stats ecdf
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
#' Filter the predictors, retaining only the top \code{n}.
#'
#' This performs a univariate test of association of each predictor (not listed
#' in \code{ignore}) with the outcome, and retains the top \code{n} variables
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
