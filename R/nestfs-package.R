#' Cross-validated (nested) forward selection
#'
#' This package provides an implementation of forward selection based on linear
#' and logistic regression which adopts cross-validation as a core component of
#' the selection procedure.
#'
#' The engine of the package is \code{\link{forward.selection}}, whose aim is
#' to select a set of variables out of those available in the dataset. The
#' selection of variables can be done according to two main different criteria:
#' by paired-test p-value or by largest decrease in validation log-likelihood.
#' A combined criteria is also available.
#'
#' The role of \code{\link{nested.forward.selection}} is to allow the
#' evaluation of the selection method by providing an unbiased estimate of the
#' performance of the selected variables on withdrawn data.
#'
#' Forward selection is an inherently slow approach, as for each variable a
#' model needs to be fitted. In our implementation, this issue is further
#' aggravated by the fact that an inner cross-validation happens at each
#' iteration, with the aim of guiding the selection towards variables that
#' have better generalization properties.
#'
#' The code is parallelized over the inner folds, thanks to the \pkg{foreach}
#' package. User time therefore depends on the number of available cores, but
#' there is no advantage in using more cores than inner folds. A parallel
#' backend must be registered before starting, otherwise operations will run
#' sequentially with a warning reported. This can be done through a call to
#' \code{\link[doParallel]{registerDoParallel}}, for example.
#'
#' The main advantage of forward selection is that it provides an immediately
#' interpretable model, and the panel of variables obtained is in some sense
#' the least redundant one, particularly if the number of variables to choose
#' from is not too large (in our experience, up to about 30-40 variables).
#'
#' However, when the number of variables is much larger than that, forward
#' selection, besides being unbearably slow, may be more subject to
#' overfitting, which is in the nature of its greedy-like design. These
#' undesirable effects can be somewhat remedied by applying some filtering
#' (see \code{"num.filter"} argument to \code{\link{forward.selection}}, thus
#' reducing the number or variables entering the selection phase.
#'
#' @author
#' Marco Colombo \email{m.colombo@@ed.ac.uk}
#'
#' @docType package
"_PACKAGE"

#' @importFrom foreach getDoParWorkers getDoParRegistered
.onAttach <- function(libname, pkgname) {
  if (getDoParRegistered())
    packageStartupMessage("nestfs: currently using ", getDoParWorkers(),
                          " cores.")
  else
    packageStartupMessage("nestfs: register a parallel backend before starting.")
}
