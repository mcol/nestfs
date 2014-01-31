\name{nested.glm}
\alias{nested.glm}

\title{
  Cross-validated Generalized Linear Models
}
\description{
  Run linear or logistic regression on the set of cross-validation folds.
}
\usage{
  nested.glm(x, y, folds, family = c("binomial", "gaussian"))
}
\arguments{
  \item{x}{
    Data frame.
  }
  \item{y}{
    Outcome variable.
  }
  \item{folds}{
    Set of cross-validation folds.
  }
  \item{family}{
    Type of model fitted.
  }
}
\details{
  This can be used to establish a baseline model, often built only on
  the initial set of covariates (those that would be passed through the
  \code{init.vars} argument to \code{forward.selection}).
}
\value{
  A list of length equal to \code{length(folds)}, where each entry
  contains the following fields:
  \item{fit}{
    Predicted values for the withdrawn observations.
  }
  \item{obs}{
    Observed values for the withdrawn observations.
  }
  \item{llk}{
    Log-likelihood.
  }
  \item{test.idx}{
    Indices of the the withdrawn observations for this fold.
  }
}

\examples{
\dontrun{
  data(diabetes)
  all.folds <- create.folds(10, nrow(diabetes), seed=1)
  lr.res <- nested.glm(diabetes[, c("age", "sex", "bmi", "tc",
                                    "ldl", "hdl", "ltg", "glu")],
                       diabetes$Y, all.folds, family="gaussian")
  }
}

\keyword{internal}
