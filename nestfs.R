## run the forward selection starting from a set of variables or a model
forward.selection <- function(x.all, y.all, init.vars, test=c("t", "wilcoxon"),
                              sel.crit=c("paired.test", "total.loglik", "both"),
                              num.filter=0, filter.ignore=init.vars,
                              num.inner.folds=50, max.iters=30, max.pval=0.5,
                              min.llk.diff=0, seed=50,
                              init.model=NULL) {
  univ.logreg <- function(model, x.train, x.test) {
    regr <- glm(as.formula(model), data=x.train, family="binomial")
    y.pred <- predict(regr, newdata=x.test, type="response")
    y.test <- x.test$y

    loglik <- sum(log(y.pred[y.test == 1])) + sum(log(1 - y.pred[y.test == 0]))

    res <- cbind(coefficients(summary(regr))[, c(1, 4)], loglik)
    colnames(res) <- c("coef", "p.value", "valid.llk")
    return(res)
  }
  inner.fold <- function(x.all, y.all, model, other.vars, test.idx) {
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    x.train <- cbind(y=y.all[train.idx], x.all[train.idx, ])
    x.test <- cbind(y=y.all[test.idx], x.all[test.idx, ])

    ## current model
    tt.curr <- univ.logreg(model, x.train, x.test)
    all.stats <- tail(tt.curr, n=1)

    ## models augmented with one additional variable at a time
    for (var in other.vars) {
      model.var <- paste(model, var, sep=" + ")
      tt <- univ.logreg(model.var, x.train, x.test)
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
    fmt <- "%2d %35s %9.5f %9.2f %7.2f\n"
    if (iter == 1)
      cat(sprintf("%2s %35s %9s %9s %7s\n",
                  "#", "Variable", "P-value", "Log-Lik", "Diff"))
    cat(sprintf(fmt, iter, var, pval, llk, diff.llk))
  }
  stopifnot(all.equal(names(table(y.all)), c("0", "1")))
  pval.test <- match.arg(test)
  sel.crit <- match.arg(sel.crit)
  if (is.null(init.model))
    init.model <- paste("y ~", paste(init.vars, collapse= " + "))
  else {
    ## work out the variables from the initialization model
    cat("Using init.model, ignoring init.vars\n")
    stopifnot(length(grep("~", init.model)) > 0)
    init.vars <- unlist(strsplit(init.model, "~" ))[2]
    init.vars <- gsub(" ", "", init.vars)
    init.vars <- unlist(strsplit(init.vars, "\\+" ))
  }

  all.folds <- create.folds(num.inner.folds, nrow(x.all), seed=seed)
  num.init.vars <- length(init.vars)
  model.vars <- init.vars
  model.llks <- c(rep(NA, num.init.vars - 1), 0)
  model.pvals <- model.iter <- rep(NA, num.init.vars)
  model <- init.model

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
                              total.llk=rowSums(all.llk[-1, ]))

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
      univ.diffs <- all.llk[-1, ]
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

  return(list(fs=data.frame(vars=model.vars, pvals=model.pvals, llks=model.llks,
                diffs=c(NA, diff(model.llks)), iter=model.iter,
                row.names=NULL, stringsAsFactors=FALSE),
              panel=setdiff(model.vars, init.vars),
              iter1=iter1, all.iter=all.iter))
}

nested.forward.selection <- function(x.all, y.all, model.vars, all.folds,
                                     ...) {
  all.res <- list()
  num.folds <- length(all.folds)
  for (fold in 1:num.folds) {

    cat("* Outer Fold", fold, "of", num.folds,
        "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- all.folds[[fold]]
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    x.train <- x.all[train.idx, ]; y.train <- y.all[train.idx]

    fs <- forward.selection(x.train, y.train, model.vars, ...)
    this.fold <- list(test.idx)
    model <- plain.logreg(x.all[, fs$fs$vars], y.all, this.fold)[[1]]
    stopifnot(all.equal(model$caseness.test, y.all[test.idx]))
    panel <- fs$panel
    fs$fs$coef <- NA
    fs$fs$coef[match(panel, fs$fs$vars)] <- model$regr$coef[panel]
    res <- list(fs=fs$fs, fit=model$fit, caseness.test=model$caseness.test,
                panel=panel, model=summary(model$regr), test.idx=test.idx,
                iter1=fs$iter1, all.iter=fs$all.iter, call=match.call())
    all.res[[fold]] <- res
  }
  class(all.res) <- "nestfs"
  return(all.res)
}

plain.logreg <- function(x, y, folds) {
  stopifnot(all.equal(nrow(x), length(y)))
  res <- list()
  for (fold in 1:length(folds)) {
    if (fold %% 10 == 0)
      cat("Fold", fold, "\n")
    idx.test <- folds[[fold]]
    idx.train <- setdiff(1:nrow(x), idx.test)
    x.test <- x[idx.test, ]; x.train <- x[idx.train, ]
    y.test <- y[idx.test];   y.train <- y[idx.train]
    model <- paste("y.train ~", paste(colnames(x.train), collapse=" + "))
    regr <- glm(as.formula(model), data=x.train, family="binomial")
    y.pred <- predict(regr, newdata=x.test, type="response")
    loglik <- sum(log(y.pred[y.test == 1])) + sum(log(1 - y.pred[y.test == 0]))
    res[[fold]] <- list(regr=regr, fit=y.pred, caseness.test=y.test,
                        llk=loglik, test.idx=idx.test)
  }
  return(res)
}

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

## contribution to the AUCs by adding one marker at a time
aucs.fs.incremental <- function(x, y, selection, init, folds, num.keep) {
  num.folds <- length(folds)
  stopifnot(all.equal(length(selection), num.folds))
  all.res <- list()
  for (fold in 1:num.folds) {

    ## cut the selection to the specified number
    vars <- union(init, head(selection[[fold]], num.keep))

    ## fit a logistic regression model on the training/test split
    test.idx <- folds[[fold]]
    res <- plain.logreg(x[, vars], y, list(test.idx))[[1]]
    res$test.idx <- test.idx
    all.res[[fold]] <- res
  }
  return(all.res)
}
