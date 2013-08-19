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
                              min.llk.diff=0, n.add=1, rep.every=25, seed=50) {
  par.univ.logreg <- function(x.train, x.test, model, met) {
    model.met <- paste(model, met, sep=" + ")
    tt <- univ.logreg(model.met, x.train, x.test)
    return(tail(tt, n=1))
  }
  stopifnot(all.equal(names(table(y.all)), c("0", "1")))
  pval.test <- match.arg(test)

  all.folds <- produce.folds(1, num.folds, nrow(x.all), seed=seed)[[1]]
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
    diff.llk <- ifelse(iter > 1, chosen.llk - max(model.llks, na.rm=TRUE), Inf)
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
  }

  return(data.frame(vars=model.vars, pvals=model.pvals, llks=model.llks,
                    diffs=c(NA, diff(model.llks)), iter=model.iter,
                    row.names=NULL, stringsAsFactors=FALSE))
}

nested.forward.selection <- function(x.all, y.all, model.vars, all.folds,
                                     test=c("t", "wilcoxon"), num.inner.folds,
                                     max.iters=50, max.pval=0.5,
                                     min.llk.diff=0, seed=50) {
  all.res <- list()
  num.folds <- length(all.folds)
  for (fold in 1:num.folds) {

    cat("* Outer Fold", fold, "of", num.folds,
        "-", format(Sys.time(), "%H:%M"), "\n")

    test.idx <- all.folds[[fold]]
    train.idx <- setdiff(seq(nrow(x.all)), test.idx)
    x.train <- x.all[train.idx, ]; y.train <- y.all[train.idx]
    fs <- forward.selection(x.train, y.train, model.vars, test=test,
                            max.iters=max.iters, num.folds=num.inner.folds,
                            max.pval=max.pval, min.llk.diff=min.llk.diff,
                            rep.every=100, seed=seed)
    this.fold <- list(test.idx)
    model <- plain.logreg(x.all[, fs$vars], y.all, this.fold)[[1]]
    stopifnot(all.equal(model$caseness.test, y.all[test.idx]))
    res <- list(fs=fs, fit=model$fit, caseness.test=model$caseness.test,
                model=model$regr)
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

  nfolds <- length(res)
  sel <- NULL
  for (fold in 1:nfolds)
    sel <- c(sel, res[[fold]]$fs$vars)
  sel <- grep("DEMOG", sel, invert=TRUE, value=TRUE)
  tsel <- tryCatch(table(comp2bio(sel)), error=function(e) table(sel))
  ttt <- data.frame(tsel, emp.pval=1-as.numeric(tsel)/nfolds)
  ttt <- ttt[order(ttt$emp.pval), ]
  rownames(ttt) <- NULL
  return(ttt)
}
