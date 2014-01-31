\name{filter.predictors}
\alias{filter.predictors}

\title{
  Filtering of predictors
}
\description{
  Filter the predictors, retaining only the top n.
}
\usage{
  filter.predictors(x, y, n, ignore = NULL)
}
\arguments{
  \item{x}{
    Design matrix.
  }
  \item{y}{
    Outcome variable.
  }
  \item{n}{
    Number of variables to retain.
}
  \item{ignore}{
    Names of variables that should be ignored by the filter.
  }
}
\details{
  This performs a univariate test of association with the outcome, and
  retains the top \code{n} variables with smallest p-value.
}
\value{
  An array of indices of variables that are retained by the filter.
}

\keyword{internal}
