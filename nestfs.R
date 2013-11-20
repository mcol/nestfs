## run the forward selection starting from a set of variables or a model
forward.selection <- function(x.all, y.all, init.vars, test=c("t", "wilcoxon"),
                              sel.crit=c("paired.test", "total.loglik"),
                              num.filter=0, filter.ignore=init.vars,
                              num.folds=50, max.iters=30, max.pval=0.5,
                              min.llk.diff=0, n.add=1, seed=50,
                              init.model=NULL, save.iter1=NULL) {
  univ.logreg <- function(model, x.train, x.test, mean.llk=FALSE) {
    regr <- glm(as.formula(model), data=x.train, family="binomial")
    rsum <- summary(regr)

    y.pred <- predict(regr, newdata=x.test, type="response")
    y.test <- x.test$y

    acc <- sum(round(y.pred) == y.test) / length(y.pred)
    loglik <- sum(log(y.pred[y.test == 1])) + sum(log(1 - y.pred[y.test == 0]))
    if (mean.llk)
      loglik <- loglik / length(y.pred)

    res <- cbind(coef(regr), coefficients(rsum)[, 4], acc, loglik)
    colnames(res) <- c("LogEstimate", "p-value", "TestAcc", "TestLogLik")
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
    for (met in other.vars) {
      model.met <- paste(model, met, sep=" + ")
      tt <- univ.logreg(model.met, x.train, x.test)
      all.stats <- rbind(all.stats, tail(tt, n=1))
    }
    rownames(all.stats) <- c("Base", other.vars)

    return(all.stats)
  }
  paired.pvals <- function(all.llk, test=c("t", "wilcoxon")) {
    test <- match.arg(test)
    test.function <- list(t=t.test, wilcoxon=wilcox.test)
    pvals <- NULL
    for (i in 2:nrow(all.llk)) {
      ttt <- test.function[[test]](all.llk[i, ], all.llk[1, ],
                                   paired=TRUE, alternative="greater")
      pvals <- c(pvals, ttt$p.value)
    }
    names(pvals) <- rownames(all.llk)[-1]
    return(pvals)
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

  all.folds <- create.folds(num.folds, nrow(x.all), seed=seed)
  num.init.vars <- length(init.vars)
  model.vars <- init.vars
  model.llks <- c(rep(NA, num.init.vars - 1), 0)
  model.pvals <- model.iter <- rep(NA, num.init.vars)
  model <- init.model
  iter1 <- list()

  ## filtering according to association with outcome
  if (num.filter > 0) {

    ## run the filter on the training part of all inner folds
    all.filt.idx <- NULL
    for (fold in 1:num.folds) {

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

  ## variable selection
  for (iter in 1:max.iters) {

    other.vars <- setdiff(all.vars, model.vars)

    ## loop over the folds
    res.inner <- (foreach(fold=1:num.folds)
                  %dopar%
                  inner.fold(x.all, y.all, model, other.vars,
                             all.folds[[fold]]))

    ## compute the loglikelihood for the initialization model
    if (iter == 1) {
      model.llks[num.init.vars] <- sum(sapply(res.inner, function(z) z[1, 4]))
      iter1 <- res.inner
    }

    ## summarise the loglikelihoods
    all.llk <- NULL
    for (fold in 1:length(res.inner))
      all.llk <- cbind(all.llk, res.inner[[fold]][, 4])

    total.llk <- rowSums(all.llk[-1, ])
    tt.pvals <- paired.pvals(all.llk, pval.test)

    ## choose the best variable according to a paired test
    if (sel.crit == "paired.test") {
      thresh.pval <- sort(tt.pvals)[n.add]
      idx.pval <- which(tt.pvals <= thresh.pval)
      chosen.pval <- sort(tt.pvals[idx.pval])
      chosen.met <- names(chosen.pval)
      chosen.llk <- total.llk[chosen.met]
    }

    ## choose the best variable according to the total loglikelihood
    else {
      thresh.llk <- sort(total.llk, decreasing=TRUE)[n.add]
      idx.llk <- which(total.llk >= thresh.llk)
      chosen.llk <- sort(total.llk[idx.llk])
      chosen.met <- names(chosen.llk)
      chosen.pval <- tt.pvals[chosen.met]
    }

    ## report iteration summary
    diff.llk <- chosen.llk - max(model.llks, na.rm=TRUE)
    print(data.frame(chosen.pval, chosen.llk, diff.llk))

    ## check for early termination
    if (max(chosen.pval) > max.pval)
      break
    if (diff.llk < min.llk.diff)
      break

    ## append the chosen variable to the existing ones
    model.vars <- c(model.vars, chosen.met)
    model.pvals <- c(model.pvals, chosen.pval)
    model.llks <- c(model.llks, chosen.llk)
    model.iter <- c(model.iter, rep(iter, length(chosen.met)))
    model <- paste(model, chosen.met, sep=" + ")
  }

  if (!is.null(save.iter1)) save(file=save.iter1, iter1)
  return(data.frame(vars=model.vars, pvals=model.pvals, llks=model.llks,
                    diffs=c(NA, diff(model.llks)), iter=model.iter,
                    row.names=NULL, stringsAsFactors=FALSE))
}

nested.forward.selection <- function(x.all, y.all, model.vars, all.folds,
                                     test=c("t", "wilcoxon"), num.inner.folds,
                                     sel.crit=c("paired.test", "total.loglik"),
                                     max.pval=0.5, min.llk.diff=0, max.iters=50,
                                     num.filter=0, filter.ignore=model.vars,
                                     seed=50) {
  all.res <- list()
  num.folds <- length(all.folds)
  for (fold in 1:num.folds) {

    cat("* Outer Fold", fold, "of", num.folds,
        "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- all.folds[[fold]]
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    x.train <- x.all[train.idx, ]; y.train <- y.all[train.idx]

    fs <- forward.selection(x.train, y.train, model.vars, test=test,
                            sel.crit=sel.crit,
                            num.filter=num.filter, filter.ignore=filter.ignore,
                            max.iters=max.iters, num.folds=num.inner.folds,
                            max.pval=max.pval, min.llk.diff=min.llk.diff,
                            seed=seed)
    this.fold <- list(test.idx)
    model <- plain.logreg(x.all[, fs$vars], y.all, this.fold)[[1]]
    stopifnot(all.equal(model$caseness.test, y.all[test.idx]))
    res <- list(fs=fs, fit=model$fit, caseness.test=model$caseness.test,
                model=summary(model$regr))
    res$test.idx <- test.idx
    res$call <- match.call()
    all.res[[fold]] <- res
  }
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
    acc <- sum(round(y.pred) == y.test) / length(y.test)
    loglik <- sum(log(y.pred[y.test == 1])) + sum(log(1 - y.pred[y.test == 0]))
    res[[fold]] <- list(regr=regr, fit=y.pred, caseness.test=y.test,
                        acc=acc, llk=loglik, test.idx=idx.test)
  }
  return(res)
}

summary.nestfs <- function(res) {
  iqr <- function(x) quantile(x, c(0.25, 0.75))
  format.iqr <- function(x) sprintf("(%.2g, %.2g)", x[1], x[2])
  nfolds <- length(res)
  sel <- NULL
  for (fold in 1:nfolds) {
    fs <- res[[fold]]$fs[, c("vars", "diffs", "iter")]
    sel <- rbind(sel, fs[!is.na(fs$iter), ]) # exlude init model
  }
  props <- tapply(sel$vars, sel$vars, length)
  ranks <- round(tapply(sel$iter, sel$vars, median))
  ranks.iqr <- tapply(sel$iter, sel$vars, function(z) format.iqr(round(iqr(z))))
  diffs <- tapply(sel$diff, sel$vars, median)
  diffs.iqr <- tapply(sel$diff, sel$vars, function(z) format.iqr(iqr(z)))
  vars <- names(props)
  fullname <- tryCatch(comp2bio(vars), error=function(e) names(props))
  ttt <- data.frame(vars, fullname, freq=props, emp.pval=1 - props / nfolds,
                    rank=ranks, rankIQR=ranks.iqr,
                    diffLogLik=sprintf("%.3f", diffs), diffLogLikIQR=diffs.iqr)
  ttt <- ttt[order(ttt$emp.pval, ttt$rank), ]
  rownames(ttt) <- NULL
  return(ttt)
}

summary.iter1 <- function(iter1) {
  iter1.pvals <- NULL
  iter1.diffs <- NULL
  for (i in 1:length(iter1)) {
    tmp <- iter1[[i]]
    iter1.pvals <- cbind(iter1.pvals, tmp$"p-value")
    iter1.diffs <- cbind(iter1.diffs, tmp$TestLogLik - tmp$TestLogLik[1])
  }
  res <- data.frame(row.names=rownames(iter1[[1]]),
                    diffLogLik=apply(iter1.diffs, 1, sum),
                    p.value=apply(iter1.pvals, 1, median))
  res <- res[-1, ] # remove baseline
  res <- res[order(res$diffLogLik, decreasing=TRUE), ]
  return(res)
}

## contribution to the AUCs by adding one marker at a time
aucs.fs.incremental <- function(x, y, selection, clin, folds, num.keep) {
  num.folds <- length(folds)
  stopifnot(all.equal(length(selection), num.folds))
  all.res <- list()
  for (fold in 1:num.folds) {

    ## cut the selection to the specified number
    pred <- union(clin, head(selection[[fold]], num.keep))

    ## fit a logistic regression model on the training/test split
    test.idx <- folds[[fold]]
    res <- plain.logreg(x[, pred], y, list(test.idx))[[1]]
    res$test.idx <- test.idx
    all.res[[fold]] <- res
  }
  return(all.res)
}
