\name{summary.nestfs}
\alias{summary.nestfs}
\alias{summary}

\title{
  Results Summary
}
\description{
  Report summary statistics from a run of nested forward selection.
}
\usage{
  \method{summary}{nestfs}(object, iter1 = FALSE, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{"nestfs"}.
  }
  \item{iter1}{
    Whether the summary should be over all variables at the first iteration.
  }
  \item{...}{
    further arguments passed to or from other methods.
  }
}
\details{
  The function summarizes the results of a run of nested forward
  selection across the outer folds.
}
\value{
  A data frame.
}

\author{
  Marco Colombo
}
