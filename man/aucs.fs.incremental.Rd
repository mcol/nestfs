\name{aucs.fs.incremental}
\alias{aucs.fs.incremental}

\title{
  Incremental Evaluation of Nested Forward Selection
}
\description{
  This allows to evaluate the incremental contribution of each variable
  chosen by nested forward selection.
}
\usage{
  aucs.fs.incremental(x, y, selection, init, folds, num.keep)
}
\arguments{
  \item{x}{
    Data frame.
  }
  \item{y}{
    Outcome variable.
  }
  \item{selection}{
    List of variables selected in each fold.
  }
  \item{init}{
    Initial set of variables.
  }
  \item{folds}{
    Set of outer folds.
  }
  \item{num.keep}{
    Number of selected variables to evaluate.
  }
}
\details{
  To evaluate the contribution of each variable, this restricts the
  selection performed on an outer fold to the number of variables
  specified by \code{num.keep}, fits a model on the training part of the
  data, and evalutates its performance on the withdrawn observations.
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
