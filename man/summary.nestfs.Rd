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

  A function of name \code{"getfullname"} to match variable names to
  full names is searched on the current workspace, and if found full
  names are included in the summary data frame.
}
\value{
  A data frame.
}

\author{
  Marco Colombo
}
