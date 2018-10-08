context("zero_vec")
library(SIMLR)

test_that("zero_vec works correctly", {
  expect_equal(zero_vec(1:10, k = 4), c(rep(0, 6), 7, 8, 9, 10))
  expect_equal(sum(0 != zero_vec(1:100, k = 39)), 39)
  expect_equal(sum(0 != zero_vec(runif(1000), k = 39)), 39)
})
