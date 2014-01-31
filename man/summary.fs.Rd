\name{summary.fs}
\alias{summary.fs}

\title{
  Results Summary
}
\description{
  Report summary statistics from a run of a single run of forward
  selection.
}
\usage{
  \method{summary}{fs}(object, \dots)
}
\arguments{
  \item{object}{
    An object of class \code{"fs"}.
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
}
\value{
  A data frame with the following columns:

  \item{vars}{
    Variables in the initial model followed by variables selected.
  }
  \item{pvals}{
    Paired test p-values when selected.
  }
  \item{llks}{
    Validation log-likelihoods.
  }
  \item{diffs}{
    Differences in validation log-likelihoods.
  }
  \item{iter}{
    Iteration when the variable was selected.
  }
}
\note{
  A function of name \code{"getfullname"} to match variable names to
  full names is searched on the current workspace, and if found full
  names are included in the summary data frame.
}
