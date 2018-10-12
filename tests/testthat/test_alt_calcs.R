context("alt_calcs")
library(SIMLR)
library(purrr)
library(Matrix)

# Check every combination of methods on the function
method_tests <- function(methods., fn) {
  results <- lapply(methods., fn)
  names(results) <- methods.
  apply(
    combn(methods., 2),
    MARGIN = 2,
    FUN = function(m)
      expect_equal(results[[m[1]]], results[[m[2]]], check.attributes = FALSE)
  )
}


test_that("alternative methods", {
  #
  # Sample some random matrix and vector
  set.seed(3737)
  W_dense <- matrix(runif(100), nrow = 10)
  W_sparse <- Matrix(0, nrow = 10, ncol = 10)
  for (c. in 1:10) {
    W_sparse[sample(10, 3), c.] <- rnorm(3)
  }
  w <- exp(rnorm(10))
  #
  # scale_rows
  scale_row_methods <- c("dense", "sparse")
  method_tests(scale_row_methods, compose(as.matrix, partial(scale_rows, W_dense, w)))
  method_tests(scale_row_methods, compose(as.matrix, partial(scale_rows, W_sparse, w)))
  #
  # scale_cols
  scale_col_methods <- c("orig", "denom", "dense", "sparse")
  method_tests(scale_col_methods, compose(as.matrix, partial(scale_cols, W_dense, w)))
  method_tests(scale_col_methods, compose(as.matrix, partial(scale_cols, W_sparse, w)))
  #
  # col_sums
  col_sum_methods <- c("apply", "colSums")
  method_tests(col_sum_methods, compose(as.matrix, partial(col_sums, W_dense)))
  method_tests(col_sum_methods, compose(as.matrix, partial(col_sums, W_sparse)))
  #
  # row_sums
  row_sum_methods <- c("apply", "rowSums")
  method_tests(row_sum_methods, compose(as.matrix, partial(row_sums, W_dense)))
  method_tests(row_sum_methods, compose(as.matrix, partial(row_sums, W_sparse)))
  #
  # multiply_rows
  multiply_rows_methods <- c("diag", "sweep")
  method_tests(multiply_rows_methods, compose(as.matrix, partial(multiply_rows, w, W_dense)))
  method_tests(multiply_rows_methods, compose(as.matrix, partial(multiply_rows, w, W_sparse)))
  #
  # multiply_cols
  multiply_cols_methods <- c("diag", "sweep")
  method_tests(multiply_cols_methods, compose(as.matrix, partial(multiply_cols, w, W_dense)))
  method_tests(multiply_cols_methods, compose(as.matrix, partial(multiply_cols, w, W_sparse)))
})
