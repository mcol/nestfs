## create a list of indices corresponding to cross-validation folds
create.folds <- function(num.folds, num.rows, seed=NULL) {
  if (!is.null(seed))
    set.seed(seed)
  ## ensure that the test indices are not repeated
  full.idx <- numeric(ceiling(num.rows / num.folds) * num.folds)
  full.idx[1:num.rows] <- sample(num.rows)
  fold.matrix <- as.data.frame(matrix(full.idx, ncol=num.folds, byrow=TRUE))
  folds <- lapply(fold.matrix, function(z) z[z > 0])
  names(folds) <- NULL
  return(folds)
}
