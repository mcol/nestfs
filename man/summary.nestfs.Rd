\name{summary.nestfs}
\alias{summary.nestfs}

\title{
  Results Summary
}
\description{
  Report summary statistics from a run of nested forward selection
  across the outer folds.
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
    Further arguments passed to or from other methods.
  }
}
\value{
  A data frame with the following columns:

  \item{vars}{
    Variables in the initial model followed by variables selected.
  }
  \item{percent}{
    Percentage of folds in which the variable was selected.
  }
  \item{coef}{
    Median coefficient for the variable.
  }
  \item{coefIQR}{
    Inter-quartile range for the variable coefficient.
  }
  \item{rank}{
    Median iteration in which the variable was selected.
  }
  \item{rankIQR}{
    Inter-quartile range for rank of the variable.
  }
  \item{diffLogLik}{
    Median difference in log-likelihoods.
  }
  \item{diffLogLikIQR}{
    Inter-quartile range for the difference in log-likelihoods.
  }
}
\note{
  A function of name \code{"getfullname"} to match variable names to
  full names is searched on the current workspace, and if found full
  names are included in the summary data frame.
}
