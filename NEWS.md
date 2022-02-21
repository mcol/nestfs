# nestfs 1.0.1 (2022-02-21)

- Remove test that will fail due changes in R-devel.


# nestfs 1.0 (2019-09-21)

- Introduce the fs() and nested.fs() functions which adopt a new interface based
  on formulas
- Change the interface of nested.glm() to align to the new formula interface.
- Change default values for the max.iters (from 15 to 10) and min.llk.diff
  (from 0 to 2) options.
- Replace the parallel computation backend from the doParallel to the
  parallel package.
- Restructure the diabetes dataset to be a single data frame.
- Update and expand the example in the README file.
- Use markdown in the package documentation.
- Include tests in the package.


# nestfs 0.9.2 (2019-05-02)

- Silence messages output by newer versions of the pROC package.


# nestfs 0.9.1 (2018-12-16)

- Change maintainer email address.


# nestfs 0.9 (2018-09-25)

- Use getfullname() if available also in summary.fs().
- Make nested.glm() accept a formula argument so that models with interaction
  terms can be specified. This also ensures that such models are fitted
  correctly in nested.forward.selection() after selection has been performed.
- Add the family field and assign a class to the object created by
  nested.glm().
- Add nested.performance() to compute the performance of cross-validated
  models as the area under the curve or the correlation coefficient.


# nestfs 0.8.6 (2018-08-13)

- Document the default selection criterion.
- Correct the check for the verbose option in nested.forward.selection().
- Fix an error occurring in nested.forward.selection() when a categorical
  variable is selected.


# nestfs 0.8.5 (2018-08-02)

- Make the univariate filter cope with non-matching names in filter.ignore.
- Parallelise the univariate filtering step.
- Add the verbose option to forward.selection().
- Return the coefficients of summary() instead of summary() itself from
  nested.glm().
- Swap family and folds in nested.glm() for consistency with other functions.
- Add tests for nested.glm().


# nestfs 0.8.4 (2018-07-04)

- Close the parallel clusters at the end of the examples.
- Vectorize the computation of differences in log-likelihoods at iteration 1.
- First version on CRAN.


# nestfs 0.8.3 (2018-07-03)

- Rewrite the examples to satisfy the CRAN upload request.
- Decrease the minimum number of inner folds to 5.


# nestfs 0.8.2 (2018-06-29)

- Use family$dev.resids() to compute log-likelihoods.
- Fix forward.selection() when there's only one variable to choose from.
- Allow to specify variable names in the choose.from argument and not only
  indices.
- Allow more freedom in how the outcome variable can be specified for
  logistic regression.
- Rename parameters x.all, y.all and all.folds to x, y, and folds.
- Merge init.vars and init.model to make formulas a first class input type.
- Rework the diabetes dataset and save it in .rda format.
- Replace the doMC package with doParallel.
- Remove automatic registration of the parallel backend when attaching the
  package to pass checks on the R-devel win-builder machine.
- Add tests for forward.selection() and nested.forward.selection().


# nestfs 0.8.1 (2018-06-21)

- Sort the indices of the test observations within each fold.
- Reorder some arguments of forward.selection() according to importance.
- Improve the argument checks in forward.selection().
- Let the family argument also be one of the family functions.
- Add tests for argument checks.


# nestfs 0.8 (2017-07-03)

- Limit the variable names in the output to the length of the field.
- Clarify that the p-value from forward selection is a false discovery rate.
- Convert documentation to roxygen2 format.


# nestfs 0.7.20160815 (2016-08-15)

- Check that the indices in the folds don't exceed the size of the dataset.
- Make the init.model option work in more cases.


# nestfs 0.7.20150729 (2015-07-29)

- Check for missing values in the predictors and in the outcome variable.
- Return only the right-hand side of the formula in final.model from
  forward.selection().


# nestfs 0.7 (2013-12-21)

- First version of the package.
