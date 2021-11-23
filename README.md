nestfs: Cross-validated (nested) forward selection
======

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/nestfs)](https://cran.r-project.org/package=nestfs)
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
fs.res <- fs(Y ~ age + sex, diabetes, family=gaussian())
summary(fs.res)
##      vars          fdr      llks     diffs iter
## 1     age           NA        NA        NA   NA
## 2     sex           NA -2136.854        NA   NA
## 3     ltg 7.008928e-10 -2058.831 78.022766    1
## 4     bmi 1.850715e-05 -2009.568 49.263387    2
## 5     map 2.020038e-03 -1999.253 10.314799    3
## 6 age.sex 1.508210e-02 -1992.544  6.709064    4
## 7     hdl 4.039276e-02 -1985.208  7.336623    5
## 8 bmi.map 7.474167e-02 -1980.913  4.294736    6
```

By default, selection happens over all variables present in the data.frame
that are not part of the initial model. This can be controlled through the
`choose.from` option, which accepts variable names or indices.

It is possible to promote sparser selection by requesting a larger improvement
in log-likelihood (option `min.llk.diff`): this is advisable especially when the
number of variables to choose from exceeds 10-15, as it's our experience that
even the default setting of 2 (already stricter than what adopted by other
packages) may lead to some overfitting. In any case, it's possible to set a
maximum size of the panel selected by reducing the number of iterations (option
`max.iters`, by default set to 10).

Nested forward selection is helpful to assess the stability of the selection
process by performing it on each training split of the cross-validation folds:

```r
folds <- create.folds(10, nrow(diabetes), seed=1)
nest.res <- nested.fs(Y ~ age + sex, diabetes, family=gaussian(), folds=folds)
summary(nest.res)
##       vars percent    coef          coefIQR rank      rankIQR diffLogLik  diffLogLikIQR
## 1      bmi     100  24.547   (23.61, 25.48)    2 (1.00, 2.00)     61.021 (44.49, 76.85)
## 2      ltg     100  23.729   (22.39, 24.41)    2 (1.00, 2.00)     52.868 (36.09, 69.36)
## 3      map     100  15.147   (14.45, 15.88)    3 (3.00, 3.75)      8.366   (8.04, 9.61)
## 4      hdl     100 -13.297 (-13.65, -12.55)    4 (4.00, 4.00)      6.728   (6.35, 7.83)
## 5  age.sex      80   8.825     (8.72, 9.24)    5 (5.00, 6.00)      4.625   (4.45, 5.37)
## 6  bmi.map      70   8.165     (7.55, 8.27)    6 (5.50, 7.00)      3.604   (2.66, 4.15)
## 7  bmi.glu      20   4.460     (4.07, 4.85)    5 (5.00, 5.00)      3.535   (3.09, 3.98)
## 8    glu.2      20   6.477     (6.47, 6.49)    6 (6.25, 6.75)      2.984   (2.56, 3.41)
## 9  sex.map      20   6.862     (6.71, 7.01)    6 (5.25, 5.75)      2.936   (2.89, 2.98)
## 10 age.glu      10   7.469     (7.47, 7.47)    3 (3.00, 3.00)      4.826   (4.83, 4.83)
## 11 age.map      10   7.365     (7.36, 7.36)    6 (6.00, 6.00)      2.679   (2.68, 2.68)
## 12   bmi.2      10   7.987     (7.99, 7.99)    6 (6.00, 6.00)      2.466   (2.47, 2.47)
```

The output above shows that `bmi`, `ltg`, `map` and `hdl` are chosen in all
folds, and most of the improvement in fit is provided by the first two variables,
which agrees with what was found when running forward selection on all data.

Most importantly, nested forward selection produces a cross-validated measure
of performance of the selection process, which is an unbiased estimate of the
predictive performance of the selected panels on withdrawn data:
```r
nested.performance(nest.res)
## Correlation coefficient: 0.7097
```

This can be compared with what is obtained by the baseline model on the same
set of cross-validation folds:

```r
base.res <- nested.glm(Y ~ age + sex, diabetes, family=gaussian(), folds=folds)
nested.performance(base.res)
## Correlation coefficient: 0.1551
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
