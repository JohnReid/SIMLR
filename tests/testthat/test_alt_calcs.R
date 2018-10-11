context("alt_calcs")
library(SIMLR)
library(purrr)
library(Matrix)

scale_tests <- function(W, scale_fn, scale_methods) {
  w <- exp(rnorm(10))
  scaled <- lapply(scale_methods, compose(as.matrix, partial(scale_fn, W, w)))
  names(scaled) <- scale_methods
  apply(
    combn(scale_methods, 2),
    MARGIN = 2,
    FUN = function(m)
      expect_equal(scaled[[m[1]]], scaled[[m[2]]], check.attributes = FALSE))
}


test_that("scale alternatives", {
  set.seed(3737)
  W_dense <- matrix(runif(100), nrow = 10)
  W_sparse <- Matrix(0, nrow = 10, ncol = 10)
  for( c. in 1:10 ) {
    W_sparse[sample(10, 3), c.] <- rnorm(3)
  }
  #
  # scale_rows
  scale_tests(W_dense, scale_rows, c('dense', 'sparse'))
  scale_tests(W_sparse, scale_rows, c('dense', 'sparse'))
  #
  # scale_cols
  scale_tests(W_dense, scale_cols, c('orig', 'denom', 'dense', 'sparse'))
  scale_tests(W_sparse, scale_cols, c('orig', 'denom', 'dense', 'sparse'))
})
