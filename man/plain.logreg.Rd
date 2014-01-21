\name{plain.logreg}
\alias{plain.logreg}

\title{
  Cross-validated Logistic Regression
}
\description{
  Run logistic regression on the set of cross-validation folds.
}
\usage{
  plain.logreg(x, y, folds, family = c("binomial", "gaussian"))
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
    Type of link function to be used in the model.
  }
}
\details{
  This can be used to establish a baseline model, often built only on
  the initial set of covariates (those that would be passed through the
  \code{init.vars} argument to \code{forward.selection}.
}
\value{
  A set of results.
}
\author{
  Marco Colombo
}

\examples{
\dontrun{
  data(diabetes)
  all.folds <- create.folds(10, nrow(diabetes), seed=1)
  lr.res <- plain.logreg(diabetes[, c("age", "sex", "bmi", "tc",
                                      "ldl", "hdl", "ltg", "glu")],
                         diabetes$Y, all.folds, family="gaussian")
  }
}

\keyword{internal}
