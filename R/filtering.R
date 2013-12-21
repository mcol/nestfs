library(dgof)

## compare p(x|y=0) and p(x|y=1) with the Kolmogorov-Smirnov test
ks.pval <- function(x, y) {
  require(dgof)
  is.discrete <- apply(x, 1, function(z) length(unique(z)) < 10)
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

## filter the predictors, retaining only the top n
filter.predictors <- function(x, y, n, ignore=NULL) {
  stopifnot(n <= ncol(x))
  if (!is.matrix(x))
    x <- sapply(x, as.numeric)
  scores <- ks.pval(t(x[, !colnames(x) %in% ignore]), y)
  limit <- sort(scores, decreasing=TRUE)[n]
  keep.idx <- match(rownames(scores)[scores >= limit], colnames(x))
  return(keep.idx)
}
