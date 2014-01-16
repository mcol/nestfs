\name{summary.nestfs}
\alias{summary.nestfs}
\alias{summary.fs}
\alias{summary}

\title{
  Results Summary
}
\description{
  Report summary statistics from a run of forward selection or nested
  forward selection.
}
\usage{
  \method{summary}{fs}(object, \dots)

  \method{summary}{nestfs}(object, iter1 = FALSE, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{"fs"} or \code{"nestfs"}.
  }
  \item{iter1}{
    Whether the summary should be over all variables at the first iteration.
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
}
\details{
  For an object of class \code{"nestfs"} the function summarizes the
  selected panels across the outer folds.
}
\value{
  A data frame.
}

\author{
  Marco Colombo
}
