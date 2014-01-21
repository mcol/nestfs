\name{diabetes}
\alias{diabetes}
\docType{data}

\title{
  Diabetes Data with Interactions
}
\description{
  The \code{diabetes} data consists of 442 patients in which the
  outcome variable \code{Y} is a quantitative measure of disease
  progression. The data includes 10 baseline measurements, in addition
  to 45 interactions and 9 quadratic terms, for a total of 64 variables
  for each patient.

  All predictors are transformed to have zero mean and unit variance.
}
\source{
  B. Efron, T. Hastie, I. Johnstone and R. Tibshirani,
  "Least angle regression (with discussion)",
  \emph{Ann. Statist.}, 32 (2), 407-499, 2004.

  \url{http://www.stanford.edu/~hastie/Papers/LARS/data64.txt}
}

\examples{
  data(diabetes, package="nestfs")
}

\keyword{datasets}
