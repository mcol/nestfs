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
                              num.folds=50, max.iters=30, max.pval=0.15) {
  par.univ.logreg <- function(x.train, x.test, model, met) {
    model.met <- paste(model, met, sep=" + ")
    tt <- univ.logreg(model.met, x.train, x.test)
    return(tail(tt, n=1))
  }
  stopifnot(all.equal(names(table(y.all)), c("0", "1")))
  pval.test <- match.arg(test)

  all.folds <- produce.folds(1, num.folds, nrow(x.all), seed=50)[[1]]
  all.vars <- colnames(x.all)
  model.pvals <- rep(NA, length(model.vars))
  model.llks <- rep(NA, length(model.vars))

  for (iter in 1:max.iters) {

    cat("*** Iteration", iter, "\n")

    model <- paste("y ~", paste(model.vars, collapse= " + "))
    other.vars <- setdiff(all.vars, model.vars)

    ## loop over the folds
    res.by.fold <- NULL
    for (fold in 1:length(all.folds)) {

      if (fold %% 25 == 0)
        cat("Fold", fold, "\n")

      test.idx <- all.folds[[fold]]
      train.idx <- setdiff(seq(nrow(x.all)), test.idx)
      x.train <- cbind(y=y.all[train.idx], x.all[train.idx, ])
      x.test <- cbind(y=y.all[test.idx], x.all[test.idx, ])

      ## models augmented with one metabolite at a time
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
    tt.pvals <- paired.pvals(all.llk)[[pval.test]]
    idx.min <- which.min(tt.pvals)
    chosen.pval <- tt.pvals[idx.min]
    chosen.met <- names(idx.min)
    chosen.llk <- sum(all.llk[chosen.met, ])
    cat(chosen.met, chosen.pval, chosen.llk, "\n")

    ## append the chosen variable to the existing ones
    model.vars <- c(model.vars, chosen.met)
    model.pvals <- c(model.pvals, chosen.pval)
    model.llks <- c(model.llks, chosen.llk)

    ## check for early termination
    if (chosen.pval > max.pval)
      break
  }

  return(data.frame(vars=model.vars, pvals=model.pvals, llks=model.llks))
}

paired.pvals <- function(all.llk) {
  wilcoxon.pvals <- t.pvals <- NULL
  for (i in 2:nrow(all.llk)) {
    tmp <- wilcox.test(all.llk[i, ], all.llk[1, ],
                       paired=TRUE, alternative="greater")
    wilcoxon.pvals <- c(wilcoxon.pvals, tmp$p.value)
    tmp <- t.test(all.llk[i, ], all.llk[1, ],
                  paired=TRUE, alternative="greater")
    t.pvals <- c(t.pvals, tmp$p.value)
  }
  names(wilcoxon.pvals) <- names(t.pvals) <- rownames(all.llk)[-1]
  return(list(wilcoxon=wilcoxon.pvals, t=t.pvals))
}

plain.logreg <- function(x, y, folds) {
  stopifnot(all.equal(nrow(x), length(y)))
  all.regr <- list()
  all.acc <- all.llk <- NULL
  for (fold in 1:length(folds)) {
    cat(fold, "\n")
    idx.test <- folds[[fold]]
    idx.train <- setdiff(1:nrow(x), idx.test)
    x.test <- x[idx.test, ]; x.train <- x[idx.train, ]
    y.test <- y[idx.test];   y.train <- y[idx.train]
    preds <- paste("~", paste(colnames(x.train), collapse=" + "))
    model <- paste(quote(y.train), preds)
    regr <- glm(as.formula(model), data=x.train, family="binomial")
    y.pred <- predict(regr, newdata=x.test, type="response")
    all.regr[[fold]] <- regr
    all.acc <- c(all.acc, sum(round(y.pred) == y.test) / length(y.test))
    loglik <- sum(log(y.pred[y.test == 1])) + sum(log(1 - y.pred[y.test == 0]))
    all.llk <- c(all.llk, loglik)
  }
  return(list(acc=all.acc, llk=all.llk, regr=all.regr))
}
