\name{create.folds}
\alias{create.folds}
\title{
  Cross-Validation Folds
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

\examples{
  all.folds <- create.folds(50, 307, 0)
}
