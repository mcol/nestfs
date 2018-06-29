#' Results summary for forward selection
#'
#' Report summary statistics from a single run of forward selection.
#'
#' @param object An object of class \code{fs}.
#' @param ... Further arguments passed to or from other methods.
#'        These are currently ignored.
#'
#' @return
#' A dataframe with the following columns:
#' \describe{
#' \item{vars:}{Variables in the initial model followed by variables selected.}
#' \item{fdr:}{False discovery rate, corresponding to the paired test p-values
#'       computed when the variable was selected.}
#' \item{llks:}{Validation log-likelihoods.}
#' \item{diffs:}{Differences in validation log-likelihoods.}
#' \item{iter:}{Iteration when the variable was selected.}
#' }
#'
#' @note
#' A function of name \code{"getfullname"} to match variable names to full
#' names is searched on the current workspace, and if found full names are
#' included in the summary dataframe.
#'
#' @export
summary.fs <- function(object, ...) {
  return(object$fs)
}

#' Results summary for nested forward selection
#'
#' Report summary statistics from a run of nested forward selection across the
#' outer folds.
#'
#' @param object An object of class \code{nestfs}.
#' @param iter1 Whether the summary should be over all variables at the first
#'        iteration: this can be interpreted as a cross-validated univariate
#'        test for association.
#' @param ... Further arguments passed to or from other methods.
#'        These are currently ignored.
#'
#' @return
#' A dataframe with the following columns:
#' \describe{
#' \item{vars:}{Variables in the initial model followed by variables selected.}
#' \item{percent:}{Percentage of folds in which the variable was selected.}
#' \item{coef:}{Median coefficient for the variable.}
#' \item{coefIQR:}{Inter-quartile range for the variable coefficient.}
#' \item{rank:}{Median iteration in which the variable was selected.}
#' \item{rankIQR:}{Inter-quartile range for rank of the variable.}
#' \item{diffLogLik:}{Median difference in log-likelihoods.}
#' \item{diffLogLikIQR:}{Inter-quartile range for the difference in
#' log-likelihoods.}
#' }
#'
#' @note
#' A function of name \code{"getfullname"} to match variable names to full
#' names is searched on the current workspace, and if found full names are
#' included in the summary dataframe.
#'
#' @importFrom stats median quantile
#' @export
summary.nestfs <- function(object, iter1=FALSE, ...) {
  iqr <- function(x) quantile(x, c(0.25, 0.75), na.rm=TRUE)
  format.iqr <- function(x, n=2) sprintf("(%.*f, %.*f)", n, x[1], n, x[2])
  num.folds <- length(object)
  if (iter1) {
    ## with filtering, outer folds may contain different variables
    all.iter1 <- NULL
    for (fold in 1:num.folds) {
      fold.iter1 <- object[[fold]]$iter1[, c("total.diff.llk", "p.value")]
      fold.iter1 <- data.frame(vars=rownames(fold.iter1), fold.iter1)
      all.iter1 <- rbind(all.iter1, fold.iter1)
    }

    ## summarise all variables at the first iteration
    all.vars <- all.iter1$vars
    med.diffs <- round(tapply(all.iter1$total.diff.llk, all.vars, median), 2)
    iqr.diffs <- tapply(all.iter1$total.diff.llk, all.vars,
                        function(z) format.iqr(iqr(z), 2))
    med.pvals <- round(tapply(all.iter1$p.value, all.vars, median), 3)
    iqr.pvals <- tapply(all.iter1$p.value, all.vars,
                        function(z) format.iqr(iqr(z), 3))
    vars <- names(med.diffs)
    fullname <- tryCatch(get("getfullname")(vars), error=function(e) NULL)
    res <- data.frame(row.names=vars,
                      med.diff.llk=med.diffs, iqr.diff.llk=iqr.diffs,
                      med.pvalue=med.pvals, iqr.pvalue=iqr.pvals)
    if (!is.null(fullname)) res <- cbind(fullname, res, stringsAsFactors=FALSE)
    res <- res[order(res$med.diff.llk, decreasing=TRUE), ]
    return(res)
  }
  sel <- NULL
  for (fold in 1:num.folds) {
    fs <- object[[fold]]$fs[, c("vars", "coef", "diffs", "iter")]
    sel <- rbind(sel, fs[!is.na(fs$iter), ]) # exclude init model
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

#' Display the contents of (nested) forward selection object
#'
#' Report summary statistics from a run of (nested) forward selection.
#'
#' @param x An object of class \code{nestfs} or \code{fs}.
#' @param ... Further arguments passed to or from other methods.
#'        These are currently ignored.
#'
#' @return
#' The function returns \code{summary} invisibly.
#'
#' @seealso
#' \code{\link{forward.selection}} and \code{\link{nested.forward.selection}}.
#' @export
print.nestfs <- function(x, ...) {
  print(summary(x))
}

#' @rdname print.nestfs
#' @export
print.fs <- function(x, ...) {
  print(summary(x))
}
