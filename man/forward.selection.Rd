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
                    num.filter = 0, filter.ignore = init.vars,
                    num.inner.folds = 30, max.iters = 15, max.pval = 0.5,
                    min.llk.diff = 0, seed = 50, init.model = NULL)
}
\arguments{
  \item{x.all}{
    Data frame.
  }
  \item{y.all}{
    Outcome variable.
  }
  \item{init.vars}{
    Initial set of variables.
  }
  \item{test}{
    Type of statistical test to use.
  }
  \item{family}{
    Type of link function to be used in the model.
  }
  \item{sel.crit}{
    Selection criterion: "paired.test" chooses the variable with best
    p-value on the paired test indicated by \code{test}; "total.loglik"
    chooses the variable that provides the largest increase in
    log-likelihood; "both" attempts to combine both previous criteria,
    choosing the variable that produces the largest increase in
    log-likelihood only among the best 5 variables ranked according to the
    paired-test p-value.
  }
  \item{num.filter}{
    Number of variables to be retained by the filter (0 to use all).
  }
  \item{filter.ignore}{
    Regular expression for variables that should be ignored by the filter
    (so that they are always retained).
  }
  \item{num.inner.folds}{
    Number of inner folds.
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
  This function is used to choose a panel of variables on the data. As
  it uses all observations in the \code{x.all} data frame, it is not
  possible to produce unbiased estimates of the predictive performance
  of the panel selected (use \code{nested.forward.selection} for that).
}
\value{
  A list containing the following fields:
  \item{fs}{
    A data frame containing the forward selection summary.
  }
  \item{panel}{
    Names of variables selected (in order).
  }
  \item{iter1}{
    Summary statistics for all variables at the first iteration.
  }
  \item{all.iter}{
    Validation log-likelihoods for all inner folds at all iterations.
  }
}
\author{
  Marco Colombo
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
