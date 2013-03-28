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

forward.selection <- function(x.all, y.all, model.vars, test=c("t", "wilcoxon"),
                              num.folds=50, max.iters=30, max.pval=0.15,
                              n.add=1, rep.every=25) {
  par.univ.logreg <- function(x.train, x.test, model, met) {
    model.met <- paste(model, met, sep=" + ")
    tt <- univ.logreg(model.met, x.train, x.test)
    return(tail(tt, n=1))
  }
  stopifnot(all.equal(names(table(y.all)), c("0", "1")))
  pval.test <- match.arg(test)

  all.folds <- produce.folds(1, num.folds, nrow(x.all), seed=50)[[1]]
  all.vars <- colnames(x.all)
  model.pvals <- model.llks <- model.iter <- rep(NA, length(model.vars))

  for (iter in 1:max.iters) {

    model <- paste("y ~", paste(model.vars, collapse= " + "))
    other.vars <- setdiff(all.vars, model.vars)

    ## loop over the folds
    res.by.fold <- NULL
    for (fold in 1:length(all.folds)) {

      if (fold %% rep.every == 0)
        cat("Fold", fold, "\n")

      test.idx <- all.folds[[fold]]
      train.idx <- setdiff(seq(nrow(x.all)), test.idx)
      x.train <- cbind(y=y.all[train.idx], x.all[train.idx, ])
      x.test <- cbind(y=y.all[test.idx], x.all[test.idx, ])

      ## models augmented with one additional variable at a time
      all.stats <- (foreach(idx=1:length(other.vars), .combine=rbind)
                    %dopar%
                    par.univ.logreg(x.train, x.test, model, other.vars[idx]))

      ## current model
      tt <- univ.logreg(model, x.train, x.test)
      all.stats <- rbind(tail(tt, n=1), all.stats)
      rownames(all.stats) <- c("Base", other.vars)

      ## store results for this fold
      res.by.fold[[fold]] <- all.stats
    }

    ## summarise the loglikelihoods
    all.llk <- NULL
    for (fold in 1:length(res.by.fold))
      all.llk <- cbind(all.llk, res.by.fold[[fold]][, 4])

    ## choose the best variable according to a paired test
    tt.pvals <- paired.pvals(all.llk, pval.test)
    thresh.pval <- sort(tt.pvals)[n.add]
    idx.pval <- which(tt.pvals <= thresh.pval)
    chosen.pval <- sort(tt.pvals[idx.pval])
    chosen.met <- names(chosen.pval)
    chosen.llk <- rowSums(all.llk[chosen.met, , drop=FALSE])
    print(data.frame(chosen.pval, chosen.llk))

    ## append the chosen variable to the existing ones
    model.vars <- c(model.vars, chosen.met)
    model.pvals <- c(model.pvals, chosen.pval)
    model.llks <- c(model.llks, chosen.llk)
    model.iter <- c(model.iter, rep(iter, length(chosen.met)))

    ## check for early termination
    if (max(chosen.pval) > max.pval)
      break
  }

  return(data.frame(vars=model.vars, pvals=model.pvals, llks=model.llks,
                    iter=model.iter,
                    row.names=NULL, stringsAsFactors=FALSE))
}

nested.forward.selection <- function(x.all, y.all, model.vars,
                                     test=c("t", "wilcoxon"), num.folds=50,
                                     num.inner.folds=50, max.iters=30,
                                     max.pval=0.3) {
  all.folds <- produce.folds(1, num.folds, nrow(x.all), seed=50)[[1]]
  all.res <- list()
  for (fold in 1:length(all.folds)) {

    cat("* Outer Fold", fold, "\n")

    test.idx <- all.folds[[fold]]
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    x.train <- x.all[train.idx, ]; y.train <- y.all[train.idx]
    fs <- forward.selection(x.train, y.train, model.vars, test=test,
                            max.iters=max.iters, num.folds=num.inner.folds,
                            max.pval=max.pval, rep.every=100)
    this.fold <- list(test.idx)
    ttt <- plain.logreg(x.all, y.all, this.fold)
    stopifnot(all.equal(unlist(ttt$ground), y.all[test.idx]))
    res <- list(fs=fs, fit=unlist(ttt$pred), caseness.test=unlist(ttt$ground))
    res$test.idx <- test.idx
    res$call <- match.call()
    all.res[[fold]] <- res
  }
  return(all.res)
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

plain.logreg <- function(x, y, folds) {
  stopifnot(all.equal(nrow(x), length(y)))
  all.regr <- all.pred <- all.ground <- list()
  all.acc <- all.llk <- NULL
  for (fold in 1:length(folds)) {
    if (fold %% 10 == 0)
      cat("Fold", fold, "\n")
    idx.test <- folds[[fold]]
    idx.train <- setdiff(1:nrow(x), idx.test)
    x.test <- x[idx.test, ]; x.train <- x[idx.train, ]
    y.test <- y[idx.test];   y.train <- y[idx.train]
    preds <- paste("~", paste(colnames(x.train), collapse=" + "))
    model <- paste(quote(y.train), preds)
    regr <- glm(as.formula(model), data=x.train, family="binomial")
    y.pred <- predict(regr, newdata=x.test, type="response")
    all.regr[[fold]] <- regr
    all.pred[[fold]] <- y.pred
    all.ground[[fold]] <- y.test
    all.acc <- c(all.acc, sum(round(y.pred) == y.test) / length(y.test))
    loglik <- sum(log(y.pred[y.test == 1])) + sum(log(1 - y.pred[y.test == 0]))
    all.llk <- c(all.llk, loglik)
  }
  return(list(acc=all.acc, llk=all.llk, regr=all.regr,
              pred=all.pred, ground=all.ground))
}
