#' Run forward selection
#'
#' Run forward selection starting from a set of variables or a model.
#'
#' At each iteration, this function runs cross-validation to choose which
#' variable enters the final panel by fitting the current model (defined in
#' \code{init.model} as a list of variables or a formula for the first
#' iteration) augmented by one of the remaining variables at a time.
#'
#' This function is used to choose a panel of variables on the data. As it uses
#' all observations in the \code{x} dataframe, it is not possible to
#' produce unbiased estimates of the predictive performance of the panel
#' selected (use \code{\link{nested.forward.selection}} for that purpose).
#'
#' In the case of a binary outcome when very large number of predictors is
#' available, it may be convenient to apply a univariate association filter.
#' If \code{num.filter} is set to a positive value, then all available
#' predictors (excluding those whose name is matched by \code{filter.ignore})
#' are tested for univariate association with the outcome, and only the first
#' \code{num.filter} enter the selection phase, while the others are filtered
#' out. This is done on the training part of all inner folds. Filtering can
#' enhance the performance of forward selection when the number of available
#' variables exceeds about 30-40.
#'
#' @template args-forward
#' @template args-outcome
#' @template args-family
#' @param choose.from Indices or variable names over which the selection should
#'        be performed. If \code{NULL} (default), all variables in \code{x}
#'        that are not in \code{init.model} are considered.
#' @param test Type of statistical paired test to use (ignored if
#'        \code{sel.crit="total.loglik"}).
#' @param sel.crit Selection criterion: \code{"paired.test"} chooses the
#'        variable with best p-value on the paired test indicated by
#'        \code{test}; \code{"total.loglik"} chooses the variable that provides
#'        the largest increase in log-likelihood; \code{"both"} attempts to
#'        combine both previous criteria, choosing the variable that produces
#'        the largest increase in log-likelihood only among the best 5
#'        variables ranked according to the paired-test p-value.
#' @param num.filter Number of variables to be retained by the univariate
#'        association filter (see \strong{Details}), which can only be enabled
#'        if \code{family=binomial()}. Variables listed in \code{init.model}
#'        are never filtered. If set to 0 (default), the filter is disabled.
#' @param filter.ignore Vector of variable names that should not be pruned by
#'        the univariate association filter so that they are always allowed to
#'        be selected (ignored if \code{num.filter=0}).
#' @param num.inner.folds Number of folds in the inner cross-validation. It
#'        must be at least 5 (default: 30).
#' @param max.iters Maximum number of iterations (default: 15).
#' @param max.pval Interrupt the selection when the best achievable p-value
#'        exceeds this threshold (default: 0.5).
#' @param min.llk.diff Interrupt the selection when the best achievable
#'        improvement in log-likelihood is smaller than this threshold
#'        (default: 0).
#' @param seed Seed of the random number generator for the inner folds.
#'
#' @return
#' An object of class \code{fs} containing the following fields:
#' \describe{
#' \item{fs:}{A dataframe containing the forward selection summary.}
#' \item{init:}{The set of variables used in the initial model.}
#' \item{panel:}{ Names of variables selected (in order).}
#' \item{init.model:}{Right-hand side of the formula corresponding to the
#'       initial model.}
#' \item{final.model:}{Right-hand side of the formula corresponding to the
#'       final model.}
#' \item{family:}{Type of model fitted.}
#' \item{params:}{List of parameters used.}
#' \item{iter1:}{Summary statistics for all variables at the first iteration.}
#' \item{all.iter:}{Validation log-likelihoods for all inner folds at all
#'       iterations.}
#' }
#'
#' @examples
#' # register a parallel cluster with two cores
#' library(doParallel)
#' registerDoParallel(2)
#'
#' data(diabetes)
#' fs.res <- forward.selection(X.diab, Y.diab, ~ age + sex, family=gaussian(),
#'                             choose.from=1:10, num.inner.folds=5, max.iters=3)
#' summary(fs.res)
#'
#' # close the parallel cluster
#' stopImplicitCluster()
#' @seealso \code{\link{nested.forward.selection}}
#' @keywords multivariate
#' @importFrom foreach foreach %dopar%
#' @importFrom stats coefficients glm predict t.test update wilcox.test
#' @importFrom utils head tail
#' @export
forward.selection <- function(x, y, init.model, family,
                              choose.from=NULL, test=c("t", "wilcoxon"),
                              sel.crit=c("paired.test", "total.loglik", "both"),
                              num.filter=0, filter.ignore=NULL,
                              num.inner.folds=30, max.iters=15, max.pval=0.5,
                              min.llk.diff=0, seed=50) {
  univ.glm <- function(model, xy.train, xy.test) {
    regr <- glm(model, data=xy.train, family=family)
    y.pred <- predict(regr, newdata=xy.test, type="response")
    y.test <- xy.test$nestfs_y_
    loglik <- loglikelihood(family, y.test, y.pred, summary(regr)$dispersion)
    res <- cbind(coefficients(summary(regr))[, c(1, 4), drop=FALSE], loglik)
    colnames(res) <- c("coef", "p.value", "valid.llk")
    return(res)
  }
  paired.pvals <- function(all.llk, test.function) {
    pvals <- NULL
    for (i in 2:nrow(all.llk)) {
      ttt <- test.function(all.llk[i, ], all.llk["Base", ],
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
  extract.call.params <- function() {
    args <- c("test", "sel.crit", "num.filter", "filter.ignore",
              "num.inner.folds", "max.iters", "max.pval", "min.llk.diff", "seed")
    args <- as.list(parent.env(environment()))[args]

    ## for arguments with multiple default values keep only the first, which
    ## is the one that is used
    lapply(args, head, n=1)
  }

  ## argument checks
  if (nrow(x) != length(y))
    stop("Mismatched dimensions.")
  y <- validate.outcome(y)
  if (is.null(choose.from))
    choose.from <- seq(ncol(x))
  else {
    if (is.integer(choose.from)) {
      if (min(choose.from) < 1 || max(choose.from) > ncol(x))
        stop("choose.from contains out of bound indices.")
    }
    else if (is.character(choose.from)) {
      choose.from <- match(choose.from, colnames(x))
      if (any(is.na(choose.from)))
        stop("choose.from contains names that cannot be matched.")
    }
    else
      stop("choose.from should be an integer or character vector.")
  }
  family <- validate.family(family)
  if (family$family == "binomial") {
    if (length(table(y)) != 2)
      stop("More than two classes in y with family=binomial().")
    if (any(y < 0 | y > 1))
      stop("Values in y must be between 0 and 1.")
  }
  pval.test <- list(t=t.test, wilcoxon=wilcox.test)[[match.arg(test)]]
  sel.crit <- match.arg(sel.crit)
  if (num.inner.folds < 5)
    stop("num.inner.folds should be at least 5.")
  if (max.iters < 1)
    stop("max.iters should be at least 1.")
  if (max.pval <= 0 || max.pval >= 1)
    stop("max.pval should be between 0 and 1.")
  if (min.llk.diff < 0)
    stop("min.llk.diff cannot be negative.")

  ## set up initial model
  init.model <- validate.init.model(init.model)
  init.vars <- setdiff(all.vars(init.model), "nestfs_y_")

  ## check that all variables exist in the dataframe of predictors
  var.match <- match(init.vars, colnames(x))
  if (any(is.na(var.match)))
    stop("'", paste(init.vars[is.na(var.match)], collapse="', '"),
         "' not present in x.")

  ## check that there is no missingness in the initial model
  if (any(is.na(x[, init.vars])))
    stop("Missing values in the variables of the initial model.")

  ## create the inner folds
  folds <- create.folds(num.inner.folds, nrow(x), seed=seed)

  ## if the model contains only the intercept term, count 1 variable
  num.init.vars <- max(length(init.vars), 1)
  model.vars <- init.vars
  if (length(model.vars) == 0)
    model.vars <- "<empty>"
  model.llks <- c(rep(NA, num.init.vars - 1), 0)
  model.pvals <- model.iter <- rep(NA, num.init.vars)
  model <- init.model

  ## limit the number of variables to choose from
  keep.vars <- union(choose.from, match(init.vars, colnames(x)))
  keep.vars <- keep.vars[!is.na(keep.vars)] # remove NAs from interaction terms
  x <- x[, keep.vars]

  ## filtering according to association with outcome
  if (num.filter > 0) {

    if (family$family != "binomial")
      stop("num.filter can only be used with family=binomial().")
    if (num.filter >= ncol(x))
      stop("num.filter cannot exceed the number of available predictors.")

    ## run the filter on the training part of all inner folds
    all.filt.idx <- NULL
    ignore.vars <- c(init.vars, filter.ignore)
    for (fold in 1:num.inner.folds) {

      train.idx <- setdiff(seq(nrow(x)), folds[[fold]])
      x.train <- x[train.idx, ]
      y.train <- y[train.idx]
      filt.idx <- filter.predictors(x.train, y.train, num.filter,
                                    ignore=ignore.vars)
      all.filt.idx <- c(all.filt.idx, filt.idx)
    }

    ## keep the union of the variables retained in the inner folds
    filt.idx <- sort(table(all.filt.idx), decreasing=TRUE)
    filt.idx <- as.integer(names(filt.idx))[1:num.filter]
    keep.idx <- union(match(ignore.vars, colnames(x)), filt.idx)
    x <- x[, keep.idx]
  }

  ## variable selection
  all.vars <- colnames(x)
  all.iter <- list()
  iter1 <- NULL

  for (iter in 1:max.iters) {

    other.vars <- setdiff(all.vars, model.vars)
    if (length(other.vars) == 0)
      break

    ## loop over the folds
    res.inner <- foreach(fold=1:num.inner.folds) %dopar% {

      test.idx <- folds[[fold]]
      train.idx <- setdiff(seq(nrow(x)), test.idx)
      xy.train <- cbind(nestfs_y_=y[train.idx], x[train.idx, ])
      xy.test <- cbind(nestfs_y_=y[test.idx], x[test.idx, ])

      ## current model
      curr <- univ.glm(model, xy.train, xy.test)
      all.stats <- tail(curr, n=1)

      ## models augmented with one additional variable at a time
      for (var in other.vars) {
        model.augm <- update(model, paste(". ~ .", var, sep=" + "))
        augm <- univ.glm(model.augm, xy.train, xy.test)
        all.stats <- rbind(all.stats, tail(augm, n=1))
      }
      rownames(all.stats) <- c("Base", other.vars)
      return(all.stats)
    }

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
      univ.diffs <- all.llk[-1, , drop=FALSE] - rep(all.llk["Base", ],
                                                    each=length(other.vars))
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
    model <- update(model, paste(". ~ .", chosen.var, sep=" + "))
  }

  res <- list(fs=data.frame(vars=model.vars, fdr=model.pvals, llks=model.llks,
                  diffs=c(NA, diff(model.llks)), iter=model.iter,
                  row.names=NULL, stringsAsFactors=FALSE),
              init=init.vars,
              panel=setdiff(model.vars, c(init.vars, "<empty>")),
              init.model=as.character(init.model)[3],
              final.model=as.character(model)[3],
              family=family$family,
              params=extract.call.params(),
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
#' of the selected panels on withdrawn data by running forward selection on
#' a predetermined set of folds.
#'
#' @template args-forward
#' @template args-outcome
#' @template args-family
#' @template args-folds
#' @param ... Arguments to \code{forward.selection}.
#'
#' @return
#' An object of class \code{nestfs} of length equal to
#' \code{length(folds)}, where each element is an object of class
#' \code{fs} containing the following additional fields:
#' \describe{
#' \item{fit:}{Predicted values for the withdrawn observations.}
#' \item{obs:}{Observed values for the withdrawn observations.}
#' \item{test.idx:}{Indices of the the withdrawn observations for this fold.}
#' \item{model:}{Summary of the model built using the selected panel.}
#' }
#'
#' @examples
#' # register a parallel cluster with two cores
#' library(doParallel)
#' registerDoParallel(2)
#'
#' data(diabetes)
#' folds <- create.folds(2, nrow(X.diab), seed=1)
#' nestfs.res <- nested.forward.selection(X.diab, Y.diab, ~ age + sex,
#'                                        gaussian(), folds, choose.from=1:10,
#'                                        num.inner.folds=5, max.iters=3)
#' summary(nestfs.res)
#'
#' # close the parallel cluster
#' stopImplicitCluster()
#' @seealso \code{\link{forward.selection}}
#' @keywords multivariate
#' @export
nested.forward.selection <- function(x, y, init.model, family, folds, ...) {

  all.res <- list()
  num.folds <- length(folds)
  for (fold in 1:num.folds) {

    cat("* Outer Fold", fold, "of", num.folds,
        "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- folds[[fold]]
    train.idx <- setdiff(seq(nrow(x)), test.idx)
    x.train <- x[train.idx, ]; y.train <- y[train.idx]

    fs <- forward.selection(x.train, y.train, init.model, family, ...)
    this.fold <- list(test.idx)
    model <- nested.glm(x[, fs$fs$vars], y, this.fold,
                        family=family)[[1]]
    stopifnot(all.equal(model$obs, y[test.idx]))
    panel <- fs$panel
    fs$fs$coef <- NA
    fs$fs$coef[match(panel, fs$fs$vars)] <- model$coef[panel]
    fs$fit <- model$fit
    fs$obs <- model$obs
    fs$test.idx <- test.idx
    fs$model <- model$summary
    all.res[[fold]] <- fs
  }
  class(all.res) <- "nestfs"
  return(all.res)
}

#' Cross-validated generalized linear models
#'
#' Run linear or logistic regression on a set of cross-validation folds.
#'
#' This can be used to establish a baseline model, often built only on the
#' initial set of covariates (those that would be passed through the
#' \code{init.model} argument to \code{forward.selection}).
#'
#' @param x Dataframe of predictors.
#' @template args-outcome
#' @template args-folds
#' @template args-family
#' @param store.glm Whether the object produced by \code{glm} should be
#'        stored (default: \code{FALSE}).
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
#' # register a parallel cluster with two cores
#' library(doParallel)
#' registerDoParallel(2)
#'
#' data(diabetes)
#' folds <- create.folds(10, nrow(X.diab), seed=1)
#' base.res <- nested.glm(X.diab[, c("age", "sex", "bmi", "tc",
#'                                   "ldl", "hdl", "ltg", "glu")],
#'                        Y.diab, folds, gaussian())
#'
#' # close the parallel cluster
#' stopImplicitCluster()
#' @importFrom stats as.formula glm predict
#' @keywords multivariate
#' @export
nested.glm <- function(x, y, folds, family, store.glm=FALSE) {
  stopifnot(all.equal(nrow(x), length(y)))
  stopifnot(max(unlist(folds)) <= nrow(x))
  y <- validate.outcome(y)
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

#' Log-likelihood function
#'
#' Compute the log-likelihood up to a constant.
#'
#' @template args-family
#' @param obs Vector of observed values.
#' @param fit Vector of predicted values.
#' @param disp Dispersion parameter (1 for logistic regression).
#'
#' @noRd
loglikelihood <- function(family, obs, fit, disp) {
  sum(family$dev.resids(obs, fit, -0.5 / disp) - log(disp) / 2)
}

#' Validate the outcome variable
#'
#' Ensure that the outcome variable has been specified correctly.
#'
#' @param y Outcome variable to test.
#'
#' @return
#' A valid outcome variable. The function throws an error if the outcome
#' variable cannot be used.
#'
#' @noRd
validate.outcome <- function(y) {
  if (any(is.na(y)))
    stop("Outcome variable contains missing values.", call.=FALSE)
  if (is.character(y))
    stop("Outcome variable cannot be a character vector.", call.=FALSE)
  if (is.factor(y)) {
    if (length(levels(y)) != 2)
      stop("A factor outcome variable can only have two levels.", call.=FALSE)
    y <- as.integer(y) - 1
  }
  if (!(is.numeric(y) || is.integer(y) || is.logical(y)))
    stop("Outcome variable of invalid type.", call.=FALSE)
  return(y)
}

#' Validate initial model
#'
#' Ensure that the initial model has been specified correctly.
#'
#' @param model Model definition to test.
#'
#' @return
#' A formula describing the initial model. The function throws an error if the
#' model parameter cannot be used.
#'
#' @importFrom methods is
#' @importFrom stats as.formula update
#' @noRd
validate.init.model <- function(model) {
  if (is.null(model) || length(model) == 0) {
    model <- y ~ 1
  }
  else if (is.character(model)) {
    if (length(model) == 1 && grepl("~", model))
      model <- as.formula(model)
    else
      model <- as.formula(paste("y ~", paste(model, collapse= " + ")))
  }
  else if (!is(model, "formula"))
    stop("init.model specified incorrectly.")

  ## rename the left-hand side or add it if not present
  model <- update(model, "nestfs_y_ ~ .")

  return(model)
}

#' Validate the family argument
#'
#' Ensure that the family argument has been specified correctly.
#' This is inspired by code in \code{\link{glm}}.
#'
#' @param family Family argument to test.
#'
#' @return
#' A valid family. The function throws an error if the family argument cannot
#' be used.
#'
#' @importFrom methods is
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
