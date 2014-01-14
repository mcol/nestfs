\name{plain.logreg}
\alias{plain.logreg}

\title{
  Cross-validated Logistic Regression
}
\description{
  Run logistic regression on the set of cross-validation folds.
}
\usage{
  plain.logreg(x, y, folds, family = c("binomial", "gaussian"))
}
\arguments{
  \item{x}{
    Data frame.
  }
  \item{y}{
    Outcome variable.
  }
  \item{folds}{
    Set of cross-validation folds.
  }
  \item{family}{
    Type of link function to be used in the model.
  }
}
\details{
  This can be used to establish a baseline model, often built only on
  the initial set of covariates (those that would be passed through the
  \var{init.vars} argument to \code{forward.selection()}.
}
\value{
  A set of results.
}
\author{
  Marco Colombo
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

\keyword{multivariate}
