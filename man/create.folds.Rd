\name{create.folds}
\alias{create.folds}
\title{
  Cross-validation Folds
}
\description{
  Create a list of indices corresponding to cross-validation folds.
}
\usage{
  create.folds(num.folds, num.rows, seed = NULL)
}
\arguments{
  \item{num.folds}{
    Number of folds to be created.
  }
  \item{num.rows}{
    Number of observations in the dataset.
  }
  \item{seed}{
    Seed of the random number generator.
  }
}
\value{
  A list of length \code{num.folds} containing the indices of the
  observations to be withdrawn for validation in each fold.

  Note that the number of observations withdrawn in each fold may not be
  exactly the same if \code{num.folds} is not an integer divisor of
  \code{num.rows}.
}
\author{
  Marco Colombo
}

\examples{
  all.folds <- create.folds(50, 307, 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
