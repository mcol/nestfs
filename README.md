nestfs: Cross-validated (nested) forward selection
======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/nestfs)](https://cran.r-project.org/package=nestfs)
[![CRAN\_Downloads\_Badge](https://cranlogs.r-pkg.org/badges/nestfs)](https://cran.r-project.org/package=nestfs)

**nestfs** provides an implementation of forward selection based on linear
and logistic regression which adopts cross-validation as a core component of
the selection procedure.

Forward selection is an inherently slow approach, as for each variable a
model needs to be fitted. In our implementation, this issue is further
aggravated by the fact that an inner cross-validation happens at each
iteration, with the aim of guiding the selection towards variables that
have better generalization properties.

The code is parallelized over the inner folds, thanks to the **parallel**
package. User time therefore depends on the number of available cores, but
there is no advantage in using more cores than inner folds. The number of
cores assigned to computations must be registerd before starting by setting
the "mc.cores" option.

The main advantage of forward selection is that it provides an immediately
interpretable model, and the panel of variables obtained is in some sense
the least redundant one, particularly if the number of variables to choose
from is not too large (in our experience, up to about 30-40 variables).

However, when the number of variables is much larger than that, forward
selection, besides being unbearably slow, may be more subject to overfitting,
which is in the nature of its greedy-like design.

A precompiled package is
[available on CRAN](https://cran.r-project.org/package=nestfs).

## Usage

First load the package and register the number of cores to use by setting the
`mc.cores` option. If you are lucky enough to work on a large multicore machine,
best performance is achieved by registering as many cores as the number of inner
folds being used (the default is 30).

```r
library(nestfs)
options(mc.cores=10)
```

To run forward selection from a baseline model that contains only age and sex,
the following is enough:

```r
data(diabetes)
X <- diabetes[, -match("Y", colnames(diabetes))]
fs.res <- forward.selection(X, diabetes$Y, ~ age + sex, family=gaussian())
summary(fs.res)
```

By default, selection happens over all variables present in the data.frame
that are not part of the initial model. This can be controlled through the
`choose.from` option.

It is possible to promote sparser selection by requesting a larger improvement
in log-likelihood (option `min.llk.diff`, by default set to 0), or reducing the
number of iterations (option `max.iters`, by default set to 15).

To obtain a cross-validated measure of performance of the selection process,
nested forward selection should be run:

```r
cv.folds <- create.folds(10, nrow(X), seed=1)
nestfs.res <- nested.forward.selection(X, diabetes$Y, ~ age + sex,
                                       family=gaussian(), folds=cv.folds)
summary(nestfs.res)
nested.performance(nestfs.res)
```

## References

* M. Colombo, H.C. Looker, B. Farran et al.,
  Serum kidney injury molecule 1 and beta-2-microglobulin perform as well as
  larger panels for prediction of rapid decline in renal function in type 2
  diabetes, Diabetologia (2019) 62 (1): 156-168.
  https://doi.org/10.1007/s00125-018-4741-9

* H.C. Looker, M. Colombo, S. Hess et al.,
  Biomarkers of rapid chronic kidney disease progression in type 2 diabetes,
  Kidney International (2015), 88 (4): 888-896.
  https://doi.org/10.1038/ki.2015.199

* H.C. Looker, M. Colombo, F. Agakov et al.,
  Protein biomarkers for the prediction of cardiovascular disease in type 2
  diabetes, Diabetologia (2015) 58 (6): 1363-1371.
  https://doi.org/10.1007/s00125-015-3535-6
