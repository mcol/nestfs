## report the forward selection panel
summary.fs <- function(object, ...) {
  return(object$fs)
}

## summarise the results across the nested folds
summary.nestfs <- function(object, iter1=FALSE, ...) {
  iqr <- function(x) quantile(x, c(0.25, 0.75), na.rm=TRUE)
  format.iqr <- function(x, n=2) sprintf("(%.*f, %.*f)", n, x[1], n, x[2])
  num.folds <- length(object)
  sel <- NULL
  iter1.pvals <- NULL
  iter1.diffs <- NULL
  for (fold in 1:num.folds) {
    fs <- object[[fold]]$fs[, c("vars", "coef", "diffs", "iter")]
    sel <- rbind(sel, fs[!is.na(fs$iter), ]) # exlude init model
    fold.iter1 <- object[[fold]]$iter1
    iter1.diffs <- cbind(iter1.diffs, fold.iter1$total.diff.llk)
    iter1.pvals <- cbind(iter1.pvals, fold.iter1$p.value)
  }
  if (iter1) {
    ## summarise all variables at the first iteration
    med.diffs <- round(apply(iter1.diffs, 1, median), 2)
    med.pvals <- round(apply(iter1.pvals, 1, median), 3)
    iqr.diffs <- apply(iter1.diffs, 1, function(z) format.iqr(iqr(z)))
    iqr.pvals <- apply(iter1.pvals, 1, function(z) format.iqr(iqr(z), 3))
    res <- data.frame(row.names=rownames(fold.iter1),
                      med.diff.llk=med.diffs, iqr.diff.llk=iqr.diffs,
                      med.pvalue=med.pvals, iqr.pvalue=iqr.pvals)
    res <- res[order(res$med.diff.llk, decreasing=TRUE), ]
    return(res)
  }
  props <- tapply(sel$vars, sel$vars, length)
  coefs <- round(tapply(sel$coef, sel$vars, median), 3)
  coefs.iqr <- tapply(sel$coef, sel$vars, function(z) format.iqr(iqr(z)))
  ranks <- round(tapply(sel$iter, sel$vars, median))
  ranks.iqr <- tapply(sel$iter, sel$vars, function(z) format.iqr(round(iqr(z))))
  diffs <- tapply(sel$diff, sel$vars, median)
  diffs.iqr <- tapply(sel$diff, sel$vars, function(z) format.iqr(iqr(z)))
  vars <- names(props)
  fullname <- tryCatch(get("getfullname")(vars), error=function(e) NULL)
  res <- data.frame(percent=round(props/num.folds * 100, 2),
                    coef=coefs, coefIQR=coefs.iqr,
                    rank=ranks, rankIQR=ranks.iqr,
                    diffLogLik=sprintf("%.3f", diffs), diffLogLikIQR=diffs.iqr)
  if (!is.null(fullname)) res <- cbind(fullname, res, stringsAsFactors=FALSE)
  res <- cbind(vars, res, stringsAsFactors=FALSE)
  res <- res[order(-res$percent, res$rank), ]
  rownames(res) <- NULL
  return(res)
}
