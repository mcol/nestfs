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


#' Results summary for forward selection
#'
#' Report summary statistics from a single run of forward selection.
#'
#' @param object,x An object of class `fs`.
#' @param ... Further arguments passed to or from other methods.
#'        These are currently ignored.
#'
#' @return
#' A data frame with the following columns:
#' \item{vars}{Variables in the initial model followed by variables selected.}
#' \item{fdr}{False discovery rate, corresponding to the paired test p-values
#'       computed when the variable was selected.}
#' \item{llks}{Validation log-likelihoods.}
#' \item{diffs}{Differences in validation log-likelihoods.}
#' \item{iter}{Iteration when the variable was selected.}
#'
#' @note
#' A function of name `"getfullname"` to match variable names to full
#' names is searched on the current workspace, and if found full names are
#' included in the summary data frame.
#'
#' @export
summary.fs <- function(object, ...) {
  res <- object$fs
  vars <- res$vars
  fullname <- tryCatch(get("getfullname")(vars), error=function(e) NULL)
  if (!is.null(fullname))
    res <- cbind(vars, fullname, res[, -1], stringsAsFactors=FALSE)
  return(res)
}

#' @rdname summary.fs
#' @export
print.fs <- function(x, ...) {
  print(summary(x))
}

#' Results summary for nested forward selection
#'
#' Report summary statistics from a run of nested forward selection across the
#' outer folds.
#'
#' @param object,x An object of class `nestfs`.
#' @param iter1 Whether the summary should be over all variables at the first
#'        iteration: this can be interpreted as a cross-validated univariate
#'        test for association.
#' @param ... Further arguments passed to or from other methods.
#'        These are currently ignored.
#'
#' @return
#' A data frame with the following columns:
#' \item{vars}{Variables selected.}
#' \item{percent}{Percentage of folds in which the variable was selected.}
#' \item{coef}{Median coefficient for the variable.}
#' \item{coefIQR}{Inter-quartile range for the variable coefficient.}
#' \item{rank}{Median iteration in which the variable was selected.}
#' \item{rankIQR}{Inter-quartile range for rank of the variable.}
#' \item{diffLogLik}{Median difference in log-likelihoods.}
#' \item{diffLogLikIQR}{Inter-quartile range for the difference in
#'       log-likelihoods.}
#'
#' @note
#' A function of name `"getfullname"` to match variable names to full
#' names is searched on the current workspace, and if found full names are
#' included in the summary data frame.
#'
#' @export
summary.nestfs <- function(object, iter1=FALSE, ...) {
  format.iqr <- function(x, n=2) {
    x <- quantile(x, c(0.25, 0.75), na.rm=TRUE)
    sprintf("(%.*f, %.*f)", n, x[1], n, x[2])
  }
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
    med.diffs <- tapply(all.iter1$total.diff.llk, all.vars, median)
    iqr.diffs <- tapply(all.iter1$total.diff.llk, all.vars, format.iqr)
    med.pvals <- tapply(all.iter1$p.value, all.vars, median)
    iqr.pvals <- tapply(all.iter1$p.value, all.vars,
                        function(z) format.iqr(z, 3))
    vars <- names(med.diffs)
    fullname <- tryCatch(get("getfullname")(vars), error=function(e) NULL)
    res <- data.frame(row.names=vars,
                      med.diff.llk=round(med.diffs, 2),
                      iqr.diff.llk=iqr.diffs,
                      med.pvalue=round(med.pvals, 3),
                      iqr.pvalue=iqr.pvals)
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
  coefs <- tapply(sel$coef, sel$vars, median)
  ranks <- tapply(sel$iter, sel$vars, median)
  diffs <- tapply(sel$diff, sel$vars, median)
  coefs.iqr <- tapply(sel$coef, sel$vars, format.iqr)
  ranks.iqr <- tapply(sel$iter, sel$vars, function(z) format.iqr(round(z)))
  diffs.iqr <- tapply(sel$diff, sel$vars, format.iqr)
  vars <- names(props)
  fullname <- tryCatch(get("getfullname")(vars), error=function(e) NULL)
  res <- data.frame(percent=round(props/num.folds * 100, 2),
                    coef=round(coefs, 3), coefIQR=coefs.iqr,
                    rank=round(ranks, 0), rankIQR=ranks.iqr,
                    diffLogLik=round(diffs, 3), diffLogLikIQR=diffs.iqr)
  if (!is.null(fullname)) res <- cbind(fullname, res, stringsAsFactors=FALSE)
  res <- cbind(vars, res, stringsAsFactors=FALSE)
  res <- res[order(-res$percent, res$rank), ]
  rownames(res) <- NULL
  return(res)
}

#' @rdname summary.nestfs
#' @export
print.nestfs <- function(x, ...) {
  print(summary(x))
}

#' Compute cross-validated performance
#'
#' Compute an unbiased estimate of the performance of a given model or
#' forward selected panel using the results obtained on the cross-validation
#' folds.
#'
#' @param x An object of class `nestfs` or `nestglm`.
#'
#' @return
#' An object of class `nestperf` containing the following fields:
#' \item{observed}{Vector of observed values from all folds.}
#' \item{predicted}{Vector of predicted values from all folds.}
#' \item{performance}{A performance measure: the area under the curve (AUC) if
#'       `family="binomial"`, or the correlation coefficient if
#'       `family="gaussian"`.}
#'
#' @seealso [nested.fs()] and [nested.glm()].
#' @importFrom pROC auc
#' @export
nested.performance <- function(x) {
  if (!inherits(x, c("nestfs", "nestglm")))
    stop("Object is not of 'nestfs' or 'nestglm' class.", call.=FALSE)

  ## summarise observed and predicted values
  num.folds <- length(x)
  obs <- fit <- NULL
  for (fold in 1:num.folds) {
    res <- x[[fold]]
    obs <- c(obs, res$obs)
    fit <- c(fit, res$fit)
  }

  ## compute auc or correlation coefficient
  is.auc <- x[[1]]$family == "binomial"
  per <- ifelse(is.auc, as.numeric(auc(obs, fit, direction="<", quiet=TRUE)),
                        cor(obs, fit))
  res <- list(observed=obs, predicted=fit, performance=per)
  attr(res, "measure") <- ifelse(is.auc, "auc", "correlation")
  class(res) <- "nestperf"
  return(res)
}

#' @rdname nested.performance
#'
#' @param digits Number of significant figures to print.
#' @param ... Further arguments passed to or from other methods.
#'        These are currently ignored.
#'
#' @export
print.nestperf <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  msg <- ifelse(attr(x, "measure") == "auc",
                "Area under the curve: ", "Correlation coefficient: ")
  cat(msg, signif(x$performance, digits=digits), "\n", sep="")
  invisible(x)
}
