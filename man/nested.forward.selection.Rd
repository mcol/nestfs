\name{nested.forward.selection}
\alias{nested.forward.selection}

\title{
  Nested Forward Selection
}
\description{
  Run nested forward selection starting from a set of variables or a model.
}
\usage{
  nested.forward.selection(x.all, y.all, init.vars, all.folds,
                           family = c("binomial", "gaussian"), ...)
}
\arguments{
  \item{x.all}{
    Data frame.
  }
  \item{y.all}{
    Outcome variable.
  }
  \item{init.vars}{
    Initial set of variables.
  }
  \item{all.folds}{
    Set of outer folds.
  }
  \item{family}{
    Type of link function to be used in the model.
  }
  \item{...}{
    Arguments to \code{forward.selection}.
  }
}
\details{
  This function allows to obtain an unbiased estimate of the performance
  of the selected panels on withdrawn data.
}
\value{
  An object of class "nestfs" of length equal to
  \code{length(all.folds)}, where each element is an object of class
  \code{"fs"} containing the following additional fields:
  \item{fit}{
    Predicted values for the withdrawn observations.
  }
  \item{obs}{
    Observed values for the withdrawn observations.
  }
  \item{test.idx}{
    Indices of the the withdrawn observations for this fold.
  }
  \item{model}{
    Summary of the model built using the selected panel.
  }
}
\author{
  Marco Colombo
}
\seealso{
  \code{\link{forward.selection}}
}

\examples{
\dontrun{
  data(diabetes)
  all.folds <- create.folds(10, nrow(diabetes), seed=1)
  nfs.res <- nested.forward.selection(diabetes[, -1], diabetes$Y,
                                      c("age", "sex"), all.folds,
                                      family="gaussian")
  summary(nfs.res)
  }
}

\keyword{multivariate}
