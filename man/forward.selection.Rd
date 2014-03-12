\name{forward.selection}
\alias{forward.selection}

\title{
  Forward Selection
}
\description{
  Run forward selection starting from a set of variables or a model.
}
\usage{
  forward.selection(x.all, y.all, init.vars, test = c("t", "wilcoxon"),
                    family = c("binomial", "gaussian"),
                    sel.crit = c("paired.test", "total.loglik", "both"),
                    choose.from = seq(ncol(x.all)),
                    num.filter = 0, filter.ignore = init.vars,
                    num.inner.folds = 30, max.iters = 15, max.pval = 0.5,
                    min.llk.diff = 0, seed = 50, init.model = NULL)
}
\arguments{
  \item{x.all}{
    Data frame of predictors: this should include all variables in the
    initial set and the variables that are allowed to enter the selected
    panel.
  }
  \item{y.all}{
    Outcome variable: if \code{family="binomial"}, it has to only have
    \code{0-1} entries.
  }
  \item{init.vars}{
    Initial set of variables.
  }
  \item{test}{
    Type of statistical paired test to use: it's ignored if
    \code{sel.crit="total.loglik"}.
  }
  \item{family}{
    Type of model fitted: if \code{"binomial"}, then the outcome
    variable is expected to only have \code{0-1} entries.
  }
  \item{sel.crit}{
    Selection criterion: \code{"paired.test"} chooses the variable with
    best p-value on the paired test indicated by \code{test};
    \code{"total.loglik"} chooses the variable that provides the largest
    increase in log-likelihood; \code{"both"} attempts to combine both
    previous criteria, choosing the variable that produces the largest
    increase in log-likelihood only among the best 5 variables ranked
    according to the paired-test p-value.
  }
  \item{choose.from}{
    Indices of the variables among which the selection should be done.
  }
  \item{num.filter}{
    Number of variables to be retained by the filter (0 to use all).
  }
  \item{filter.ignore}{
    Regular expression for variables that should be ignored by the filter
    (so that they are always retained).
  }
  \item{num.inner.folds}{
    Number of folds in the inner cross-validation.
  }
  \item{max.iters}{
    Maximum number of iterations.
  }
  \item{max.pval}{
    Interrupt the selection when the best achievable p-value exceeds
    this threshold.
  }
  \item{min.llk.diff}{
    Interrupt the selection when the best achievable improvement in
    log-likelihood is smaller than this threshold.
  }
  \item{seed}{
    Seed of the random number generator.
  }
  \item{init.model}{
    Formula that describes the initial model: if specified, this
    overrides \code{init.vars}.
  }
}
\details{
  At each iteration, this function runs cross-validation to choose which
  variable enters the final panel by fitting the current model (defined
  by \code{init.vars} at the first iteration) augmented by one of the
  remaining variables at a time.

  This function is used to choose a panel of variables on the data. As
  it uses all observations in the \code{x.all} data frame, it is not
  possible to produce unbiased estimates of the predictive performance
  of the panel selected (use \code{\link{nested.forward.selection}} for
  that purpose).

  If \code{num.filter} is positive, then all available predictors
  (excluding those whose name is matched by \code{filter.ignore}) are
  tested for univariate association with the outcome. This is done on
  the training part of all inner folds, and the top \code{num.filter}
  are retained for selection, while the others are filtered
  out. Filtering can enhance the performance of forward selection when
  the number of available variables exceeds about 30--40.
}
\value{
  An object of class \code{"fs"} containing the following fields:
  \item{fs}{
    A data frame containing the forward selection summary.
  }
  \item{init}{
    The set of variables used in the initialization.
  }
  \item{panel}{
    Names of variables selected (in order).
  }
  \item{family}{
    Type of model fitted.
  }
  \item{call}{
    The call that created this object.
  }
  \item{iter1}{
    Summary statistics for all variables at the first iteration.
  }
  \item{all.iter}{
    Validation log-likelihoods for all inner folds at all iterations.
  }
}
\seealso{
  \code{\link{nested.forward.selection}}
}

\examples{
\dontrun{
  data(diabetes)
  fs.res <- forward.selection(diabetes[, -1], diabetes$Y,
                              c("age", "sex"), family="gaussian",
                              max.iters=5)
  summary(fs.res)
  }
}

\keyword{multivariate}
