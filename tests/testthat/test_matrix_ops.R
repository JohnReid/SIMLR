context("matrix_ops")
library(SIMLR)

makesym <- function(X) (X + t(X)) / 2

test_that("diffusion operator", {
  set.seed(37)
  N <- 171
  dense <- matrix(rnorm(N^2), nrow=N)
  dense_sy <- makesym(dense)
  Dinv <- Diagonal(x = 1 / colSums(dense_sy))
  expect_equal(dn(dense_sy, 'ave'), Dinv %*% dense_sy, check.attributes = FALSE)
})
