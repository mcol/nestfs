#' Run forward selection
#'
#' Run forward selection starting from a set of variables or a model.
#'
#' At each iteration, this function runs cross-validation to choose which
#' variable enters the final panel by fitting the current model (defined either
#' by \code{init.vars} or \code{init.model} at the first iteration) augmented
#' by one of the remaining variables at a time.
#'
#' This function is used to choose a panel of variables on the data. As it uses
#' all observations in the \code{x.all} dataframe, it is not possible to
#' produce unbiased estimates of the predictive performance of the panel
#' selected (use \code{\link{nested.forward.selection}} for that purpose).
#'
#' If \code{num.filter} is positive, then all available predictors (excluding
#' those whose name is matched by \code{filter.ignore}) are tested for
#' univariate association with the outcome. This is done on the training part
#' of all inner folds, and the top \code{num.filter} are retained for
#' selection, while the others are filtered out. Filtering can enhance the
#' performance of forward selection when the number of available variables
#' exceeds about 30--40.
#'
#' @template args-forward
#' @template args-family
#' @param choose.from Indices of the variables among which the selection should
#'        be done.
#' @param test Type of statistical paired test to use (ignored if
#'        \code{sel.crit="total.loglik"}).
#' @param sel.crit Selection criterion: \code{"paired.test"} chooses the
#'        variable with best p-value on the paired test indicated by
#'        \code{test}; \code{"total.loglik"} chooses the variable that provides
#'        the largest increase in log-likelihood; \code{"both"} attempts to
#'        combine both previous criteria, choosing the variable that produces
#'        the largest increase in log-likelihood only among the best 5
#'        variables ranked according to the paired-test p-value.
#' @param num.filter Number of variables to be retained by the filter (0 to use
#'        all).
#' @param filter.ignore Regular expression for variables that should be ignored
#'        by the filter (so that they are always retained).
#' @param num.inner.folds Number of folds in the inner cross-validation.
#' @param max.iters Maximum number of iterations.
#' @param max.pval Interrupt the selection when the best achievable p-value
#'        exceeds this threshold.
#' @param min.llk.diff Interrupt the selection when the best achievable
#'        improvement in log-likelihood is smaller than this threshold.
#' @param seed Seed of the random number generator for the inner folds.
#' @param init.model Formula that describes the initial model, where the
#'        outcome variable should be called \code{y}: if specified, this
#'        overrides \code{init.vars}.
#'
#' @return
#' An object of class \code{"fs"} containing the following fields:
#' \describe{
#' \item{fs:}{A dataframe containing the forward selection summary.}
#' \item{init:}{The set of variables used in the initialization.}
#' \item{panel:}{ Names of variables selected (in order).}
#' \item{final.model:}{Right-hand side of the formula corresponding to the
#'       final model.}
#' \item{family:}{Type of model fitted.}
#' \item{call:}{The call that created this object.}
#' \item{iter1:}{Summary statistics for all variables at the first iteration.}
#' \item{all.iter:}{Validation log-likelihoods for all inner folds at all
#'       iterations.}
#' }
#'
#' @examples
#' \dontrun{
#' data(diabetes)
#' fs.res <- forward.selection(diabetes[, -1], diabetes$Y,
#'                             c("age", "sex"), family=gaussian(),
#'                             max.iters=5)
#' summary(fs.res)
#'
#' # using a formula for the initial model
#' fs.res.0 <- forward.selection(diabetes[, -1], diabetes$Y,
#'                               init.model="y ~ 1", family=gaussian(),
#'                               max.iters=5)
#' }
#' @seealso \code{\link{nested.forward.selection}}
#' @keywords multivariate
#' @importFrom foreach foreach %dopar%
#' @export
forward.selection <- function(x.all, y.all, init.vars, family,
                              choose.from=seq(ncol(x.all)), test=c("t", "wilcoxon"),
                              sel.crit=c("paired.test", "total.loglik", "both"),
                              num.filter=0, filter.ignore=init.vars,
                              num.inner.folds=30, max.iters=15, max.pval=0.5,
                              min.llk.diff=0, seed=50,
                              init.model=NULL) {
  univ.glm <- function(model, xy.train, xy.test) {
    regr <- glm(as.formula(model), data=xy.train, family=family)
    y.pred <- predict(regr, newdata=xy.test, type="response")
    y.test <- xy.test$y
    loglik <- loglikelihood(family, y.test, y.pred, summary(regr)$dispersion)
    res <- cbind(coefficients(summary(regr))[, c(1, 4), drop=FALSE], loglik)
    colnames(res) <- c("coef", "p.value", "valid.llk")
    return(res)
  }
  inner.fold <- function(x.all, y.all, model, other.vars, test.idx) {
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    xy.train <- cbind(y=y.all[train.idx], x.all[train.idx, ])
    xy.test <- cbind(y=y.all[test.idx], x.all[test.idx, ])

    ## current model
    tt.curr <- univ.glm(model, xy.train, xy.test)
    all.stats <- tail(tt.curr, n=1)

    ## models augmented with one additional variable at a time
    for (var in other.vars) {
      model.var <- paste(model, var, sep=" + ")
      tt <- univ.glm(model.var, xy.train, xy.test)
      all.stats <- rbind(all.stats, tail(tt, n=1))
    }
    rownames(all.stats) <- c("Base", other.vars)
    return(all.stats)
  }
  paired.pvals <- function(all.llk, test=c("t", "wilcoxon")) {
    test.function <- list(t=t.test, wilcoxon=wilcox.test)
    pvals <- NULL
    for (i in 2:nrow(all.llk)) {
      ttt <- test.function[[test]](all.llk[i, ], all.llk["Base", ],
                                   paired=TRUE, alternative="greater")
      pvals <- c(pvals, ttt$p.value)
    }
    names(pvals) <- rownames(all.llk)[-1]
    return(pvals)
  }
  report.iter <- function(iter, var, pval, llk, diff.llk) {
    fmt <- "%2d %40s %9.5f %9.2f %7.2f\n"
    if (iter == 1)
      cat(sprintf("%2s %40s %9s %9s %7s\n",
                  "#", "Variable", "FDR", "Log-Lik", "Diff"))
    cat(sprintf(fmt, iter, substr(var, 1, 40), pval, llk, diff.llk))
  }

  ## argument checks
  if (nrow(x.all) != length(y.all))
    stop("Mismatched dimensions.")
  if (any(is.na(y.all)))
    stop("Outcome variable contains missing values.")
  if (min(choose.from) < 1 || max(choose.from) > ncol(x.all))
    stop("choose.from contains out of bound indices.")
  family <- validate.family(family)
  if (family$family == "binomial")
    stopifnot(all.equal(names(table(y.all)), c("0", "1")))
  pval.test <- match.arg(test)
  sel.crit <- match.arg(sel.crit)
  if (num.inner.folds < 10)
    stop("num.inner.folds should be at least 10.")
  if (max.iters < 1)
    stop("max.iters should be at least 1.")
  if (max.pval <= 0 || max.pval >= 1)
    stop("max.pval should be between 0 and 1.")
  if (min.llk.diff < 0)
    stop("min.llk.diff cannot be negative.")

  if (is.null(init.model)) {
    stopifnot(all(init.vars %in% colnames(x.all)))
    init.model <- paste("y ~", paste(init.vars, collapse= " + "))
  }
  else {
    if (!missing(init.vars))
      cat("Using init.model, ignoring init.vars\n")

    ## work out the variables from the initialization model
    stopifnot(length(grep("~", init.model)) > 0)
    model.terms <- terms(as.formula(init.model))
    init.vars <- attributes(model.terms)$term.labels
  }

  ## check that there is no missingness in the variables of the initial model,
  ## excluding the interaction terms
  stopifnot(all(!is.na(x.all[, init.vars[!grepl(":", init.vars)]])))

  ## create the inner folds
  all.folds <- create.folds(num.inner.folds, nrow(x.all), seed=seed)

  ## if the model contains only the intercept term, count 1 variable
  num.init.vars <- max(length(init.vars), 1)
  model.vars <- init.vars
  if (length(model.vars) == 0)
    model.vars <- "<empty>"
  model.llks <- c(rep(NA, num.init.vars - 1), 0)
  model.pvals <- model.iter <- rep(NA, num.init.vars)
  model <- init.model

  ## limit the number of variables to choose from
  keep.vars <- union(choose.from, match(init.vars, colnames(x.all)))
  keep.vars <- keep.vars[!is.na(keep.vars)] # remove NAs from interaction terms
  x.all <- x.all[, keep.vars]

  ## filtering according to association with outcome
  if (num.filter > 0) {

    ## run the filter on the training part of all inner folds
    all.filt.idx <- NULL
    for (fold in 1:num.inner.folds) {

      train.idx <- setdiff(seq(nrow(x.all)), all.folds[[fold]])
      x.train <- x.all[train.idx, ]
      y.train <- y.all[train.idx]
      filt.idx <- filter.predictors(x.train, y.train, num.filter,
                                    ignore=filter.ignore)
      all.filt.idx <- c(all.filt.idx, filt.idx)
    }

    ## keep the union of the variables retained in the inner folds
    filt.idx <- sort(table(all.filt.idx), decreasing=TRUE)
    filt.idx <- as.integer(names(filt.idx))[1:num.filter]
    keep.idx <- union(match(filter.ignore, colnames(x.train)), filt.idx)
    x.all <- x.all[, keep.idx]
  }
  all.vars <- colnames(x.all)
  all.iter <- list()

  ## variable selection
  for (iter in 1:max.iters) {

    other.vars <- setdiff(all.vars, model.vars)
    if (length(other.vars) == 0)
      break

    ## loop over the folds
    res.inner <- (foreach(fold=1:num.inner.folds)
                  %dopar%
                  inner.fold(x.all, y.all, model, other.vars,
                             all.folds[[fold]]))

    ## collect all validation log-likelihoods
    all.llk <- NULL
    for (fold in 1:num.inner.folds)
      all.llk <- cbind(all.llk, res.inner[[fold]][, 3])
    all.iter[[iter]] <- all.llk
    inner.stats <- data.frame(p.value=paired.pvals(all.llk, pval.test),
                              total.llk=rowSums(all.llk[-1, , drop=FALSE]))

    ## choose the best variable according to a paired test
    if (sel.crit == "paired.test")
      idx.sel <- which.min(inner.stats$p.value)

    ## choose the best variable according to the total log-likelihood
    if (sel.crit == "total.loglik")
      idx.sel <- which.max(inner.stats$total.llk)

    ## choose the best variables according to the total log-likelihood among
    ## those with smallest paired p-value
    else if (sel.crit == "both") {
      ord <- order(inner.stats$p.value)[1:5]
      idx.sel <- ord[which.max(inner.stats$total.llk[ord])]
    }

    ## apply selection
    chosen.var <- rownames(inner.stats)[idx.sel]
    chosen.llk <- inner.stats$total.llk[idx.sel]
    chosen.pval <- inner.stats$p.value[idx.sel]
    if (iter == 1) {
      ## compute the log-likelihood for the initialization model
      model.llks[num.init.vars] <- sum(all.llk["Base", ])

      ## differences in validation log-likelihoods
      univ.diffs <- all.llk[-1, , drop=FALSE]
      for (i in 1:nrow(univ.diffs))
        univ.diffs[i, ] <- univ.diffs[i, ] - all.llk["Base", ]
      iter1 <- data.frame(median.diff.llk=apply(univ.diffs, 1, median),
                          total.diff.llk=apply(univ.diffs, 1, sum),
                          p.value=inner.stats$p.value)
    }

    ## report iteration summary
    diff.llk <- chosen.llk - max(model.llks, na.rm=TRUE)
    report.iter(iter, chosen.var, chosen.pval, chosen.llk, diff.llk)

    ## check for early termination
    if (max(chosen.pval) > max.pval)
      break
    if (diff.llk < min.llk.diff)
      break

    ## append the chosen variable to the existing ones
    model.vars <- c(model.vars, chosen.var)
    model.pvals <- c(model.pvals, chosen.pval)
    model.llks <- c(model.llks, chosen.llk)
    model.iter <- c(model.iter, rep(iter, length(chosen.var)))
    model <- paste(model, chosen.var, sep=" + ")
  }

  res <- list(fs=data.frame(vars=model.vars, fdr=model.pvals, llks=model.llks,
                  diffs=c(NA, diff(model.llks)), iter=model.iter,
                  row.names=NULL, stringsAsFactors=FALSE),
              init=init.vars,
              panel=setdiff(model.vars, init.vars),
              final.model=gsub("^y ~ ", "", model),
              family=family$family,
              call=match.call(),
              iter1=iter1,
              all.iter=all.iter)
  class(res) <- "fs"
  return(res)
}

#' Run nested forward selection
#'
#' Run nested forward selection starting from a set of variables or a model.
#'
#' This function allows to obtain an unbiased estimate of the performance
#' of the selected panels on withdrawn data.
#'
#' @template args-forward
#' @param all.folds Set of cross-validation folds.
#' @param ... Arguments to \code{forward.selection}.
#'
#' @return
#' An object of class \code{"nestfs"} of length equal to
#' \code{length(all.folds)}, where each element is an object of class
#' \code{"fs"} containing the following additional fields:
#' \describe{
#' \item{fit:}{Predicted values for the withdrawn observations.}
#' \item{obs:}{Observed values for the withdrawn observations.}
#' \item{test.idx:}{Indices of the the withdrawn observations for this fold.}
#' \item{model:}{Summary of the model built using the selected panel.}
#' }
#'
#' @examples
#' \dontrun{
#' data(diabetes)
#' all.folds <- create.folds(10, nrow(diabetes), seed=1)
#' nestfs.res <- nested.forward.selection(diabetes[, -1], diabetes$Y,
#'                                        c("age", "sex"), all.folds,
#'                                        family=gaussian())
#' summary(nestfs.res)
#' }
#' @seealso \code{\link{forward.selection}}
#' @keywords multivariate
#' @export
nested.forward.selection <- function(x.all, y.all, init.vars, all.folds, ...) {

  family <- list(...)$family
  all.res <- list()
  num.folds <- length(all.folds)
  for (fold in 1:num.folds) {

    cat("* Outer Fold", fold, "of", num.folds,
        "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- all.folds[[fold]]
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    x.train <- x.all[train.idx, ]; y.train <- y.all[train.idx]

    fs <- forward.selection(x.train, y.train, init.vars, ...)
    this.fold <- list(test.idx)
    model <- nested.glm(x.all[, fs$fs$vars], y.all, this.fold,
                        family=family)[[1]]
    stopifnot(all.equal(model$obs, y.all[test.idx]))
    panel <- fs$panel
    fs$fs$coef <- NA
    fs$fs$coef[match(panel, fs$fs$vars)] <- model$coef[panel]
    fs$fit <- model$fit
    fs$obs <- model$obs
    fs$test.idx <- test.idx
    fs$model <- model$summary
    fs$call <- match.call()
    all.res[[fold]] <- fs
  }
  class(all.res) <- "nestfs"
  return(all.res)
}

#' Cross-validated generalized linear models
#'
#' Run linear or logistic regression on a set of cross-validation folds
#'
#' This can be used to establish a baseline model, often built only on the
#' initial set of covariates (those that would be passed through the
#' \code{init.vars} argument to \code{forward.selection}).
#'
#' @param x Dataframe of predictors.
#' @param y Outcome variable.
#' @param folds Set of cross-validation folds.
#' @template args-family
#' @param store.glm Whether the object produced by \code{glm} should be
#'        stored.
#'
#' @return
#' A list of length equal to \code{length(folds)}, where each entry contains
#' the following fields:
#' \describe{
#' \item{summary:}{Summary of the fitted model.}
#' \item{coef:}{Coefficients of the fitted model.}
#' \item{fit:}{Predicted values for the withdrawn observations.}
#' \item{obs:}{Observed values for the withdrawn observations.}
#' \item{test.llk:}{Test log-likelihood.}
#' \item{test.idx:}{Indices of the the withdrawn observations for this fold.}
#' \item{regr:}{Object created by glm (only if \code{store.glm=TRUE}).}
#' }
#'
#' @examples
#' \dontrun{
#' data(diabetes)
#' all.folds <- create.folds(10, nrow(diabetes), seed=1)
#' base.res <- nested.glm(diabetes[, c("age", "sex", "bmi", "tc",
#'                                     "ldl", "hdl", "ltg", "glu")],
#'                        diabetes$Y, all.folds, family=gaussian())
#' }
#' @export
nested.glm <- function(x, y, folds, family, store.glm=FALSE) {
  stopifnot(all.equal(nrow(x), length(y)))
  stopifnot(max(unlist(folds)) <= nrow(x))
  family <- validate.family(family)
  res <- list()
  for (fold in 1:length(folds)) {
    if (fold %% 10 == 0)
      cat("Fold", fold, "\n")
    idx.test <- folds[[fold]]
    idx.train <- setdiff(1:nrow(x), idx.test)
    x.test <- x[idx.test, ]; x.train <- x[idx.train, ]
    y.test <- y[idx.test];   y.train <- y[idx.train]
    model <- paste("y.train ~", paste(colnames(x.train), collapse=" + "))
    regr <- glm(as.formula(model), data=x.train, family=family)
    y.pred <- predict(regr, newdata=x.test, type="response")
    loglik <- loglikelihood(family, y.test, y.pred, summary(regr)$dispersion)
    res[[fold]] <- list(summary=summary(regr), coef=regr$coef,
                        fit=y.pred, obs=y.test,
                        test.llk=loglik, test.idx=idx.test)
    if (store.glm) res[[fold]]$regr <- regr
  }
  return(res)
}

#' Log-likelihood function.
#'
#' Compute the log-likelihood up to a constant.
#'
#' @param family Type of model fitted.
#' @param obs Vector of observed values.
#' @param fit Vector of predicted values.
#' @param disp Dispersion parameter (1 for logistic regression).
#'
#' @keywords internal
loglikelihood <- function(family, obs, fit, disp) {
  sum(family$dev.resids(obs, fit, -0.5 / disp) - log(disp) / 2)
}

#' This is inspired by code in \code{\link{glm}}.
#' @noRd
validate.family <- function(family) {
  if (missing(family))
    stop("Argument of 'family' is missing.", call.=FALSE)
  if (is.character(family))
    tryCatch(
      family <- get(family, mode="function", envir=parent.frame(2)),
      error=function(e)
        stop("'", family, "' is not a valid family.", call.=FALSE)
    )
  if (is.function(family))
    family <- family()
  if (!is(family, "family"))
    stop("Argument of 'family' is not a valid family.", call.=FALSE)
  if (!family$family %in% c("gaussian", "binomial"))
    stop("Only supported families are 'gaussian' or 'binomial'.", call.=FALSE)

  return(family)
}
