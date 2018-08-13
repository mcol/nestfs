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
#' By default variables are selected according to the \code{paired.test}
#' criterion. At each iteration, the sampling distribution of differences in
#' validation log-likelihood obtained across all inner cross-validation folds
#' of the models with and without each additional variable are tested against
#' the null hypothesis of zero mean (with the alternative hypothesis being
#' that the model with the additional variable is better). The test is paired
#' according to the inner folds. Although the training folds are not
#' independent, the p-value from this test approximates the probability that
#' including the marker will not decrease the validation log-likelihood
#' (approximate false discovery rate).
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
#'        variable with smallest p-value using the paired test specified by
#'        \code{test} (see \strong{Details}), as long as this is smaller than
#'        \code{max.pval}; \code{"total.loglik"} picks the variable that gives
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
#' @param verbose Whether the variable chosen at each iteration should be
#'        printed out (default: \code{TRUE}).
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
                              min.llk.diff=0, seed=50, verbose=TRUE) {
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
    return(args)
  }

  ## argument checks
  if (nrow(x) != length(y))
    stop("Mismatched dimensions.")
  y <- validate.outcome(y)
  choose.from <- validate.choose.from(choose.from, x)
  family <- validate.family(family, y)
  test <- match.arg(test)
  pval.test <- list(t=t.test, wilcoxon=wilcox.test)[[test]]
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
  if (anyNA(var.match))
    stop("'", paste(init.vars[is.na(var.match)], collapse="', '"),
         "' not present in x.")

  ## check that there is no missingness in the initial model
  if (anyNA(x[, init.vars]))
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
  x <- x[, keep.vars]

  ## filtering according to association with outcome
  if (num.filter > 0) {

    if (family$family != "binomial")
      stop("num.filter can only be used with family=binomial().")
    if (num.filter >= ncol(x))
      stop("num.filter cannot exceed the number of available predictors.")
    if (!is.null(filter.ignore) && !is.character(filter.ignore))
      stop("filter.ignore should be a character vector or NULL.")

    ## run the filter on the training part of all inner folds
    fold <- NULL   # silence a note raised by R CMD check
    ignore.vars <- c(init.vars, filter.ignore)
    all.filt.idx <- foreach(fold=1:num.inner.folds, .combine=c) %dopar% {

      train.idx <- setdiff(seq(nrow(x)), folds[[fold]])
      x.train <- x[train.idx, ]
      y.train <- y[train.idx]
      filt.idx <- filter.predictors(x.train, y.train, num.filter,
                                    ignore=ignore.vars)
    }

    ## keep the union of the variables retained in the inner folds
    filt.idx <- sort(table(all.filt.idx), decreasing=TRUE)
    filt.idx <- as.integer(names(filt.idx))[1:num.filter]
    keep.idx <- union(match(ignore.vars, colnames(x)), filt.idx)
    keep.idx <- keep.idx[!is.na(keep.idx)]
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
    all.llk <- sapply(res.inner, function(z) z[, "valid.llk"])
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
    if (verbose)
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

  res <- list()
  family <- validate.family(family, y)
  folds <- validate.folds(folds, x)
  verbose <- isTRUE(list(...)$verbose) || is.null(list(...)$verbose)

  num.folds <- length(folds)
  for (fold in 1:num.folds) {

    if (verbose)
      cat("* Outer Fold", fold, "of", num.folds,
          "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- folds[[fold]]
    train.idx <- setdiff(seq(nrow(x)), test.idx)
    x.train <- x[train.idx, ]; y.train <- y[train.idx]

    fs <- forward.selection(x.train, y.train, init.model, family, ...)
    model <- glm.inner(x[, fs$fs$vars], y, test.idx, family)
    stopifnot(all.equal(model$obs, y[test.idx]))

    ## extract the model coefficients and report them in the forward selection
    ## object.
    ## note that factor variables will have names in model$summary that now
    ## include the level labels: by using pmatch(), we attempt a partial match
    ## for variable names. this works well if the factor had only two levels;
    ## however, if the factor had more than two levels, names cannot be matched
    ## even partially and their coefficient is set to NA (this is better than
    ## reporting for the whole variable a coefficient that corresponds only to
    ## one of the available levels).
    panel <- fs$panel
    fs$fs$coef <- NA
    idx.coefs <- pmatch(panel, rownames(model$summary))
    fs$fs$coef[match(panel, fs$fs$vars)] <- model$summary[idx.coefs, "Estimate"]
    fs$fit <- model$fit
    fs$obs <- model$obs
    fs$test.idx <- test.idx
    fs$model <- model$summary
    res[[fold]] <- fs
  }
  class(res) <- "nestfs"
  return(res)
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
#' @template args-family
#' @template args-folds
#' @param store.glm Whether the object produced by \code{glm} should be
#'        stored (default: \code{FALSE}).
#'
#' @return
#' A list of length equal to \code{length(folds)}, where each entry contains
#' the following fields:
#' \describe{
#' \item{summary:}{Summary of the coefficients of the model fitted on the
#'       training observations.}
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
#'                        Y.diab, gaussian(), folds)
#'
#' # close the parallel cluster
#' stopImplicitCluster()
#' @keywords multivariate
#' @export
nested.glm <- function(x, y, family, folds, store.glm=FALSE) {

  ## argument checks
  if (nrow(x) != length(y))
    stop("Mismatched dimensions.")
  y <- validate.outcome(y)
  family <- validate.family(family, y)
  folds <- validate.folds(folds, x)

  res <- lapply(folds, function(z) glm.inner(x, y, z, family, store.glm))
  return(res)
}


#' Fit a linear or logistic regression model on a given cross-validation fold
#'
#' Fit a model using all predictors provided on the training observations, then
#' test it on the withdrawn observations.
#'
#' @param x Dataframe of predictors containing all and only the variables to be
#'        used in the model.
#' @template args-outcome
#' @param idx.test Indices of observations to withdraw.
#' @template args-family
#' @param store.glm Whether the object produced by \code{glm} should be
#'        stored (default: \code{FALSE}).
#'
#' @return
#' A list of length equal to \code{length(folds)}.
#'
#' @importFrom stats as.formula glm predict
#' @noRd
glm.inner <- function(x, y, idx.test, family, store.glm=FALSE) {
  idx.train <- setdiff(1:nrow(x), idx.test)
  model <- paste("y ~", paste(colnames(x), collapse=" + "))
  regr <- glm(as.formula(model), data=x, family=family, subset=idx.train)
  y.pred <- predict(regr, newdata=x[idx.test, ], type="response")
  y.test <- y[idx.test]
  loglik <- loglikelihood(family, y.test, y.pred, summary(regr)$dispersion)
  res <- list(summary=coefficients(summary(regr)),
              fit=y.pred, obs=y.test, test.llk=loglik, test.idx=idx.test)
  if (store.glm) res$regr <- regr
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
  if (anyNA(y))
    stop("Outcome variable contains missing values.", call.=FALSE)
  if (is.character(y))
    stop("Outcome variable cannot be a character vector.", call.=FALSE)
  if (is.factor(y)) {
    if (length(levels(y)) != 2)
      stop("A factor outcome variable can only have two levels.", call.=FALSE)
    y <- as.integer(y) - 1
  }
  if (!(is.numeric(y) || is.logical(y)))
    stop("Outcome variable of invalid type.", call.=FALSE)

  return(as.numeric(y))
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
    if (any(model == ""))
      stop("init.model contains an empty string.", call.=FALSE)
    if (length(model) == 1 && grepl("~", model))
      model <- as.formula(model)
    else
      model <- as.formula(paste("y ~", paste(model, collapse= " + ")))
  }
  else if (!is(model, "formula"))
    stop("init.model specified incorrectly.", call.=FALSE)

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
#' @param y Outcome variable.
#'
#' @return
#' A valid family. The function throws an error if the family argument cannot
#' be used.
#'
#' @importFrom methods is
#' @noRd
validate.family <- function(family, y) {
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
    stop("Only 'gaussian' and 'binomial' are supported families.", call.=FALSE)

  if (family$family == "binomial") {
    if (length(table(y)) != 2)
      stop("y must contain two classes with family=binomial().", call.=FALSE)
    if (any(y < 0 | y > 1))
      stop("y must contain 0-1 values with family=binomial().", call.=FALSE)
  }

  return(family)
}

#' Validate the choose.from argument
#'
#' Ensure that the \code{choose.from} argument has been specified correctly.
#'
#' @param choose.from Argument to test.
#' @param x Dataframe of predictors.
#'
#' @return
#' A valid vector of variable indices. The function throws an error if the
#' argument cannot be used.
#'
#' @noRd
validate.choose.from <- function(choose.from, x) {
  if (is.null(choose.from))
    choose.from <- seq(ncol(x))
  else {
    if (is.numeric(choose.from)) {
      if (anyNA(choose.from))
        stop("choose.from contains missing values.", call.=FALSE)
      if (length(choose.from) > 0 &&
          (min(choose.from) < 1 || max(choose.from) > ncol(x)))
        stop("choose.from contains out of bounds indices.", call.=FALSE)
      if (any(choose.from != as.integer(choose.from)))
        stop("choose.from contains floating point values.", call.=FALSE)
    }
    else if (is.character(choose.from)) {
      choose.from <- match(choose.from, colnames(x))
      if (anyNA(choose.from))
        stop("choose.from contains names that cannot be matched.", call.=FALSE)
    }
    else
      stop("choose.from should be an integer or character vector.", call.=FALSE)
  }
  return(choose.from)
}

#' Validate the folds argument
#'
#' Ensure that the \code{folds} argument has been specified correctly.
#'
#' @param folds Argument to test.
#'
#' @return
#' A valid list of folds. The function throws an error if the argument cannot
#' be used.
#'
#' @noRd
validate.folds <- function(folds, x) {
  if (!is.list(folds))
    stop("folds expected to be a list.", call.=FALSE)
  all.idx <- unlist(folds)
  if (anyNA(all.idx))
    stop("folds contains missing values.", call.=FALSE)
  if (!is.numeric(all.idx))
    stop("folds contains non-numerical values.", call.=FALSE)
  if (any(all.idx != as.integer(all.idx)))
    stop("folds contains non-integer values", call.=FALSE)
  if (any(table(all.idx) > 1))
    stop("folds contains repeated indices.", call.=FALSE)
  if (any(all.idx <= 0 | all.idx > nrow(x)))
    stop("folds contains out of bounds indices.", call.=FALSE)
  return(folds)
}
