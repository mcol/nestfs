##=============================================================================
##
## Copyright (c) 2013-2019 Marco Colombo
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##=============================================================================


#' Cross-validated (nested) forward selection
#'
#' This package provides an implementation of forward selection based on linear
#' and logistic regression which adopts cross-validation as a core component of
#' the selection procedure.
#'
#' The engine of the package is [fs()], whose aim is
#' to select a set of variables out of those available in the dataset. The
#' selection of variables can be done according to two main different criteria:
#' by paired-test p-value or by largest decrease in validation log-likelihood.
#' A combined criteria is also available.
#'
#' The role of [nested.fs()] is to allow the
#' evaluation of the selection method by providing an unbiased estimate of the
#' performance of the selected variables on withdrawn data.
#'
#' Forward selection is an inherently slow approach, as for each variable a
#' model needs to be fitted. In our implementation, this issue is further
#' aggravated by the fact that an inner cross-validation happens at each
#' iteration, with the aim of guiding the selection towards variables that
#' have better generalization properties.
#'
#' The code is parallelized over the inner folds, thanks to the **parallel**
#' package. User time therefore depends on the number of available cores, but
#' there is no advantage in using more cores than inner folds. The number of
#' cores assigned to computations must be registered before starting by setting
#' the `"mc.cores"` option.
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
#' (see `"num.filter"` argument to [fs()], thus
#' reducing the number or variables entering the selection phase.
#'
#' @author
#' Marco Colombo \email{mar.colombo13@@gmail.com}
#'
#' @import stats
#' @docType package
"_PACKAGE"

.onAttach <- function(libname, pkgname) {

  ## number of cores used by default if not already set
  if (is.null(options()$mc.cores))
    options(mc.cores=min(ceiling(parallel::detectCores() / 2), 10)) # nocov

  packageStartupMessage("nestfs ", utils::packageVersion("nestfs"),
                        ": using ", options("mc.cores"),
                        " cores, set 'options(mc.cores)' to change.")
}
