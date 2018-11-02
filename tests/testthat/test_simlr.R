#
# The names we expect to see in the results
expected_names <- c('y', 'S', 'ydata', 'alphaK', 'execution.time')


#
# SIMLR
context("SIMLR")
test_that("normal result names are compliant", {
  suppressWarnings(normal <- SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0))
  expect_true(all(expected_names %in% names(normal)))
})

test_that("imputed result names are compliant", {
  suppressWarnings(if.impute <- SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, if.impute = TRUE))
  expect_true(all(expected_names %in% names(if.impute)))
})

test_that("normalized result names are compliant", {
  suppressWarnings(normalise <- SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, normalize = TRUE))
  expect_true(all(expected_names %in% names(normalise)))
})


#
# CIMLR
context("CIMLR")
test_that("result names are compliant", {
  suppressWarnings(cimlr <- CIMLR(X = GliomasReduced$in_X, c = 3, cores.ratio = 0))
  expect_true(all(expected_names %in% names(cimlr)))
})


#
# Ranking
context("SIMLR ranking")
test_that("result names are compliant", {
  suppressWarnings(ranks <- SIMLR_Feature_Ranking(A = BuettnerFlorian$results$S, X = BuettnerFlorian$in_X))
  expect_equal(names(ranks), c("pval", "aggR"))
})
