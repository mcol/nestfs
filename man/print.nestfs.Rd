\name{print.nestfs}
\alias{print.nestfs}
\alias{print.fs}

\title{
  Print an Object
}
\description{
  Report summary statistics from a run of forward selection or nested
  forward selection.
}
\usage{
  \method{print}{fs}(x, ...)

  \method{print}{nestfs}(x, ...)
}
\arguments{
  \item{x}{
    An object of class \code{"fs"} or \code{"nestfs"}.
  }
  \item{...}{
    Further arguments passed to or from other methods.
  }
}
\value{
  This returns \code{summary} invisibly.
}
\seealso{
  \code{\link{forward.selection}},
  \code{\link{nested.forward.selection}}
}
