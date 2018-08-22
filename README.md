nestfs: an R package for cross-validated (nested) forward selection.
======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/nestfs)](https://cran.r-project.org/package=nestfs)

This package provides an implementation of forward selection based on linear
and logistic regression which adopts cross-validation as a core component of
the selection procedure.

Forward selection is an inherently slow approach, as for each variable a
model needs to be fitted. In our implementation, this issue is further
aggravated by the fact that an inner cross-validation happens at each
iteration, with the aim of guiding the selection towards variables that
have better generalization properties.

The code is parallelized over the inner folds, thanks to the `foreach`
package. User time therefore depends on the number of available cores, but
there is no advantage in using more cores than inner folds. A parallel
backend must be registered before starting, otherwise operations will run
sequentially with a warning reported. This can be done through a call to
`registerDoParallel()`, for example.

The main advantage of forward selection is that it provides an immediately
interpretable model, and the panel of variables obtained is in some sense
the least redundant one, particularly if the number of variables to choose
from is not too large (in our experience, up to about 30-40 variables).

However, when the number of variables is much larger than that, forward
selection, besides being unbearably slow, may be more subject to overfitting,
which is in the nature of its greedy-like design.

A precompiled package is
[available on CRAN](https://CRAN.R-project.org/package=nestfs).
