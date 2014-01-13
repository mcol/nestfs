summary.nestfs <- function(res) {
  iqr <- function(x) quantile(x, c(0.25, 0.75), na.rm=TRUE)
  format.iqr <- function(x) sprintf("(%.2f, %.2f)", x[1], x[2])
  nfolds <- length(res)
  sel <- NULL
  for (fold in 1:nfolds) {
    fs <- res[[fold]]$fs[, c("vars", "coef", "diffs", "iter")]
    sel <- rbind(sel, fs[!is.na(fs$iter), ]) # exlude init model
  }
  props <- tapply(sel$vars, sel$vars, length)
  coefs <- round(tapply(sel$coef, sel$vars, median), 3)
  coefs.iqr <- tapply(sel$coef, sel$vars, function(z) format.iqr(iqr(z)))
  ranks <- round(tapply(sel$iter, sel$vars, median))
  ranks.iqr <- tapply(sel$iter, sel$vars, function(z) format.iqr(round(iqr(z))))
  diffs <- tapply(sel$diff, sel$vars, median)
  diffs.iqr <- tapply(sel$diff, sel$vars, function(z) format.iqr(iqr(z)))
  vars <- names(props)
  fullname <- tryCatch(comp2bio(vars), error=function(e) names(props))
  ttt <- data.frame(vars, fullname, percent=round(props/nfolds * 100, 2),
                    coef=coefs, coefIQR=coefs.iqr,
                    rank=ranks, rankIQR=ranks.iqr,
                    diffLogLik=sprintf("%.3f", diffs), diffLogLikIQR=diffs.iqr,
                    stringsAsFactors=FALSE)
  ttt <- ttt[order(-ttt$percent, ttt$rank), ]
  rownames(ttt) <- NULL
  return(ttt)
}

## summarise the first iteration information across the nested folds
summary.iter1 <- function(res) {
  iqr <- function(x) quantile(x, c(0.25, 0.75))
  format.iqr <- function(x, n) sprintf("(%.*f, %.*f)", n, x[1], n, x[2])
  stopifnot(class(res) == "nestfs")
  iter1.pvals <- NULL
  iter1.diffs <- NULL
  num.folds <- length(res)
  for (fold in 1:num.folds) {
    fold.iter1 <- res[[fold]]$iter1
    iter1.diffs <- cbind(iter1.diffs, fold.iter1$total.diff.llk)
    iter1.pvals <- cbind(iter1.pvals, fold.iter1$p.value)
  }
  med.diffs <- round(apply(iter1.diffs, 1, median), 2)
  med.pvals <- round(apply(iter1.pvals, 1, median), 3)
  iqr.diffs <- apply(iter1.diffs, 1, function(z) format.iqr(iqr(z), 2))
  iqr.pvals <- apply(iter1.pvals, 1, function(z) format.iqr(iqr(z), 3))
  res <- data.frame(row.names=rownames(fold.iter1),
                    med.diff.llk=med.diffs, iqr.diff.llk=iqr.diffs,
                    med.pvalue=med.pvals, iqr.pvalue=iqr.pvals)
  res <- res[order(res$med.diff.llk, decreasing=TRUE), ]
  return(res)
}
