% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{nested.performance}
\alias{nested.performance}
\alias{print.nestperf}
\title{Compute cross-validated performance}
\usage{
nested.performance(x)

\method{print}{nestperf}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{x}{An object of class \code{nestfs} or \code{nestglm}.}

\item{digits}{Number of significant figures to print.}

\item{...}{Further arguments passed to or from other methods.
These are currently ignored.}
}
\value{
An object of class \code{nestperf} containing the following fields:
\item{observed}{Vector of observed values from all folds.}
\item{predicted}{Vector of predicted values from all folds.}
\item{performance}{A performance measure: the area under the curve (AUC) if
\code{family="binomial"}, or the correlation coefficient if
\code{family="gaussian"}.}
}
\description{
Compute an unbiased estimate of the performance of a given model or
forward selected panel using the results obtained on the cross-validation
folds.
}
\seealso{
\code{\link[=nested.fs]{nested.fs()}} and \code{\link[=nested.glm]{nested.glm()}}.
}
