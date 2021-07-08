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


#' Cross-validated forward selection
#'
#' Run forward selection starting from a baseline model. As it uses
#' all observations in the input data frame, it is not possible to
#' produce unbiased estimates of the predictive performance of the panel
#' selected (use [nested.fs()] for that purpose).
#'
#' At each iteration, this function runs cross-validation to choose which
#' variable enters the final panel by fitting the current model augmented by
#' each remaining variable considered one at a time.
#'
#' By default variables are selected according to the `paired.test`
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
#' In the case of a binary outcome when very large number of predictors is
#' available, it may be convenient to apply a univariate association filter.
#' If `num.filter` is set to a positive value, then all available
#' predictors (excluding those whose name is matched by `filter.ignore`)
#' are tested for univariate association with the outcome, and only the first
#' `num.filter` enter the selection phase, while the others are filtered
#' out. This is done on the training part of all inner folds. Filtering can
#' enhance the performance of forward selection when the number of available
#' variables exceeds about 30-40.
#'
#' @template args-interface
#' @template args-family
#' @param choose.from Indices or variable names over which the selection should
#'        be performed. If `NULL` (default), all variables in `x` that are not
#'        in `init.model` are considered.
#' @param test Type of statistical paired test to use (ignored if
#'        `sel.crit="total.loglik"`).
#' @param sel.crit Selection criterion: `"paired.test"` chooses the
#'        variable with smallest p-value using the paired test specified by
#'        `test` (see **Details**), as long as this is smaller than
#'        `max.pval`; `"total.loglik"` picks the variable that gives
#'        the largest increase in log-likelihood; `"both"` attempts to
#'        combine both previous criteria, choosing the variable that produces
#'        the largest increase in log-likelihood only among the best 5
#'        variables ranked according to the paired-test p-value.
#' @param num.filter Number of variables to be retained by the univariate
#'        association filter (see **Details**), which can only be enabled
#'        if `family=binomial()`. Variables listed in `init.model`
#'        are never filtered. If set to 0 (default), the filter is disabled.
#' @param filter.ignore Vector of variable names that should not be pruned by
#'        the univariate association filter so that they are always allowed to
#'        be selected (ignored if `num.filter=0`).
#' @param num.inner.folds Number of folds in the inner cross-validation. It
#'        must be at least 5 (default: 30).
#' @param max.iters Maximum number of iterations (default: 10).
#' @param min.llk.diff Minimum improvement in log-likelihood required before
#'        selection is terminated (default: 2).
#' @param max.pval Interrupt the selection when the best achievable p-value
#'        exceeds this threshold (default: 0.5).
#' @param seed Seed of the random number generator for the inner folds.
#' @param verbose Whether the variable chosen at each iteration should be
#'        printed out (default: `TRUE`).
#'
#' @return
#' An object of class `fs` containing the following fields:
#' \item{fs}{A data frame containing the forward selection summary.}
#' \item{init}{The set of variables used in the initial model.}
#' \item{panel}{Names of variables selected (in order).}
#' \item{init.model}{Right-hand side of the formula corresponding to the
#'       initial model.}
#' \item{final.model}{Right-hand side of the formula corresponding to the
#'       final model after forward selection.}
#' \item{family}{Type of model fitted.}
#' \item{params}{List of parameters used.}
#' \item{iter1}{Summary statistics for all variables at the first iteration.}
#' \item{all.iter}{Validation log-likelihoods for all inner folds at all
#'       iterations.}
#'
#' @examples
#' \dontshow{oldopts <- options(mc.cores=2)}
#' data(diabetes)
#' fs.res <- fs(Y ~ age + sex, data=diabetes, family=gaussian(),
#'              choose.from=1:10, num.inner.folds=5, max.iters=3)
#' summary(fs.res)
#' \dontshow{options(oldopts)}
#'
#' @rdname forward.selection
#' @seealso [nested.fs()] and [summary.fs()].
#' @keywords multivariate
#' @export
fs <- function(formula, data, family, choose.from=NULL, test=c("t", "wilcoxon"),
               num.inner.folds=30, max.iters=10, min.llk.diff=2, max.pval=0.5,
               sel.crit=c("paired.test", "total.loglik", "both"),
               num.filter=0, filter.ignore=NULL, seed=50, verbose=TRUE) {
  univ.glm <- function(formula, data, idx.test, family) {
    res <- glm.inner(formula, data, idx.test, family, store.glm=TRUE)
    res <- cbind(coef(summary(res$regr))[, c(1, 4), drop=FALSE], res$test.llk)
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
  y <- validate.model.outcome(formula, data)
  x <- as.data.frame(data)
  init.model <- as.formula(formula)
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

  ## check that there is no missingness in the initial model
  init.vars <- all.vars(init.model)
  if (anyNA(x[, init.vars]))
    stop("Missing values in the variables of the initial model.")

  ## remove variables in the initial model from choose.from
  choose.from <- validate.choose.from(choose.from, x)
  choose.from <- setdiff(choose.from, match(init.vars, colnames(x)))

  ## remove outcome from the initial set of variables
  init.vars <- setdiff(init.vars, as.character(init.model)[2])

  ## create the inner folds
  folds <- create.folds(num.inner.folds, nrow(x), seed=seed)

  ## if the model contains only the intercept term, count 1 variable
  num.init.vars <- max(length(init.vars), 1)
  model.vars <- init.vars
  if (length(model.vars) == 0)
    model.vars <- "<empty>"
  model.llks <- c(rep(NA, num.init.vars - 1), 0)
  model.pvals <- model.iter <- rep(NA, num.init.vars)

  ## expand the interaction terms
  init.model <- update(init.model, . ~ .)
  model <- init.model

  ## limit the number of variables to choose from
  keep.vars <- union(choose.from, match(init.vars, colnames(x)))
  x <- x[, keep.vars]

  ## initialize the parallel cluster
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makePSOCKcluster(getOption("mc.cores", 1))
    on.exit(parallel::stopCluster(cl))
  }

  ## filtering according to association with outcome
  if (num.filter > 0) {

    if (family$family != "binomial")
      stop("num.filter can only be used with family=binomial().")
    if (num.filter >= ncol(x))
      stop("num.filter cannot exceed the number of available predictors.")
    if (!is.null(filter.ignore) && !is.character(filter.ignore))
      stop("filter.ignore should be a character vector or NULL.")

    ## run the filter on the training part of all inner folds
    ignore.vars <- c(init.vars, filter.ignore)
    par.fun <- function(fold) {
      train.idx <- setdiff(seq(nrow(x)), folds[[fold]])
      x.train <- x[train.idx, ]
      y.train <- y[train.idx]
      filt.idx <- filter.predictors(x.train, y.train, num.filter,
                                    ignore=ignore.vars)
    }

    if (.Platform$OS.type != "windows") {
      all.filt.idx <- parallel::mclapply(X=1:num.inner.folds, FUN=par.fun,
                                         mc.preschedule=FALSE)
    } else { # windows
      all.filt.idx <- parallel::parLapply(X=1:num.inner.folds, cl=cl,
                                          fun=par.fun)
    }

    ## keep the union of the variables retained in the inner folds
    filt.idx <- sort(table(unlist(all.filt.idx)), decreasing=TRUE)
    filt.idx <- as.integer(names(filt.idx))[1:num.filter]
    keep.idx <- union(match(ignore.vars, colnames(x)), filt.idx)
    keep.idx <- keep.idx[!is.na(keep.idx)]
    x <- x[, keep.idx]
  }

  ## variable selection
  all.vars <- setdiff(colnames(x), as.character(init.model)[2])
  all.iter <- list()
  iter1 <- NULL

  for (iter in 1:max.iters) {

    other.vars <- setdiff(all.vars, model.vars)
    if (length(other.vars) == 0)
      break

    ## loop over the folds
    par.fun <- function(fold) {
      test.idx <- folds[[fold]]

      ## current model
      curr <- univ.glm(model, data, test.idx, family)
      all.stats <- utils::tail(curr, n=1)

      ## models augmented with one additional variable at a time
      for (var in other.vars) {
        model.augm <- update(model, reformulate(c(".", var)))
        augm <- univ.glm(model.augm, data, test.idx, family)
        all.stats <- rbind(all.stats, utils::tail(augm, n=1))
      }
      rownames(all.stats) <- c("Base", other.vars)
      return(all.stats)
    }

    if (.Platform$OS.type != "windows") {
      res.inner <- parallel::mclapply(X=1:num.inner.folds, FUN=par.fun,
                                      mc.preschedule=FALSE)
    } else { # windows
      res.inner <- parallel::parLapply(X=1:num.inner.folds, cl=cl,
                                       fun=par.fun)
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
    model <- update(model, reformulate(c(".", chosen.var)))
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

#' @template args-forward
#' @template args-outcome
#' @param ... Further arguments to `fs`.
#'
#' @details
#' `forward.selection` provides the legacy interface used up to version 0.9.2.
#' It is considered discontinued, and in the future it will be deprecated and
#' eventually removed.
#'
#' @export
forward.selection <- function(x, y, init.model, family, ...) {
  if (nrow(x) != length(y))
    stop("Mismatched dimensions.")
  init.model <- validate.init.model(init.model)

  ## check that all variables exist in the data frame of predictors
  init.vars <- setdiff(all.vars(init.model), "nestfs_y_")
  var.match <- match(init.vars, colnames(x))
  if (anyNA(var.match))
    stop("'", paste(init.vars[is.na(var.match)], collapse="', '"),
         "' not present in x.")

  fs(init.model, cbind(x, nestfs_y_=y), family, ...)
}

#' Nested cross-validated forward selection
#'
#' Run nested forward selection starting from a set of variables or a model.
#'
#' This function allows to obtain an unbiased estimate of the performance
#' of the selected panels on withdrawn data by running forward selection on
#' a predetermined set of folds.
#'
#' @template args-interface
#' @template args-family
#' @template args-folds
#' @param ... Arguments to [fs()].
#'
#' @return
#' An object of class `nestfs` of length equal to `length(folds)`, where each
#' element is an object of class `fs` containing the following additional fields:
#' \item{fit}{Predicted values for the withdrawn observations.}
#' \item{obs}{Observed values for the withdrawn observations.}
#' \item{test.idx}{Indices of the the withdrawn observations for this fold.}
#' \item{model}{Summary of the model built using the selected panel.}
#'
#' @examples
#' \dontshow{oldopts <- options(mc.cores=2)}
#' data(diabetes)
#' folds <- create.folds(2, nrow(diabetes), seed=1)
#' nestfs.res <- nested.fs(Y ~ age + sex, diabetes, gaussian(), folds,
#'                         choose.from=1:10, num.inner.folds=5, max.iters=3)
#' summary(nestfs.res)
#' \dontshow{options(oldopts)}
#'
#' @rdname nested.forward.selection
#' @seealso
#' [fs()], [summary.nestfs()] and [nested.performance()].
#' @keywords multivariate
#' @export
nested.fs <- function(formula, data, family, folds, ...) {
  y <- validate.model.outcome(formula, data)
  x <- as.data.frame(data)
  family <- validate.family(family, y)
  folds <- validate.folds(folds, x)
  verbose <- isTRUE(list(...)$verbose) || is.null(list(...)$verbose)

  res <- list()
  num.folds <- length(folds)
  for (fold in 1:num.folds) {

    if (verbose)
      cat("* Outer Fold", fold, "of", num.folds,
          "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- folds[[fold]]
    fs <- fs(formula, x[-test.idx, ], family, ...)

    ## fit the final model on the training set and evaluate it on the test set
    model <- update(formula, reformulate(c(".", fs$final.model)))
    model <- glm.inner(model, x, test.idx, family)

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
    fs$obs <- unname(model$obs)
    fs$test.idx <- test.idx
    fs$model <- model$summary
    res[[fold]] <- fs
  }
  class(res) <- "nestfs"
  return(res)
}

#' @template args-forward
#' @template args-outcome
#'
#' @details
#' `nested.forward.selection` provides the legacy interface used up to version
#' 0.9.2. It is considered discontinued, and in the future it will be deprecated
#' and eventually removed.
#'
#' @export
nested.forward.selection <- function(x, y, init.model, family, folds, ...) {
  y <- validate.outcome(y)
  family <- validate.family(family, y)
  init.model <- validate.init.model(init.model)
  nested.fs(init.model, cbind(x, nestfs_y_=y), family, folds, ...)
}

#' Cross-validated generalized linear models
#'
#' Run linear or logistic regression on a set of cross-validation folds.
#' This can be used to establish a baseline model, often built only on the
#' initial set of covariates.
#'
#' @template args-interface
#' @template args-family
#' @template args-folds
#' @param store.glm Whether the object produced by `glm` should be
#'        stored (default: `FALSE`).
#'
#' @return
#' An object of class `nestglm` of length equal to `length(folds)`,
#' where each entry contains the following fields:
#' \item{summary}{Summary of the coefficients of the model fitted on the
#'       training observations.}
#' \item{family}{Type of model fitted.}
#' \item{fit}{Predicted values for the withdrawn observations.}
#' \item{obs}{Observed values for the withdrawn observations.}
#' \item{test.llk}{Test log-likelihood.}
#' \item{test.idx}{Indices of the the withdrawn observations for this fold.}
#' \item{regr}{Object created by `glm` (only if `store.glm=TRUE`).}
#'
#' @examples
#' \dontshow{oldopts <- options(mc.cores=2)}
#' data(diabetes)
#' folds <- create.folds(10, nrow(diabetes), seed=1)
#' res <- nested.glm(Y ~ age + sex + bmi + map, diabetes, gaussian(), folds)
#' \dontshow{options(oldopts)}
#'
#' @seealso [nested.performance()].
#' @keywords multivariate
#' @export
nested.glm <- function(formula, data, family, folds, store.glm=FALSE) {

  ## argument checks
  y <- validate.model.outcome(formula, data)
  x <- as.data.frame(data)
  model <- as.formula(formula)
  family <- validate.family(family, y)
  folds <- validate.folds(folds, x)

  ## check that all variables exist in the data frame of predictors
  vars <- setdiff(all.vars(model), "nestfs_y_")
  var.match <- match(vars, colnames(x))
  if (anyNA(var.match))
    stop("'", paste(vars[is.na(var.match)], collapse="', '"),
         "' not present in x.")

  res <- lapply(folds, function(z) glm.inner(model, x, z, family, store.glm))
  class(res) <- c("nestglm")
  return(res)
}


#' Fit a linear or logistic regression model on a given cross-validation fold
#'
#' Fit a model using all predictors provided on the training observations, then
#' test it on the withdrawn observations.
#'
#' @template args-interface
#' @template args-outcome
#' @param idx.test Indices of observations to withdraw.
#' @template args-family
#' @param store.glm Whether the object produced by `glm` should be
#'        stored (default: `FALSE`).
#'
#' @return
#' A list with the results of fitting a model on a specific fold.
#'
#' @noRd
glm.inner <- function(formula, data, idx.test, family, store.glm=FALSE) {
  regr <- glm(formula, data=data[-idx.test, ], family=family)
  y.pred <- predict(regr, newdata=data[idx.test, ], type="response")
  y.test <- validate.model.outcome(formula, data[idx.test, ])
  loglik <- loglikelihood(family, y.test, y.pred, summary(regr)$dispersion)
  res <- list(summary=coefficients(summary(regr)), family=family$family,
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
