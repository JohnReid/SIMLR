#
# The names we expect to see in the results
expected_names <- c('y', 'S', 'ydata', 'alphaK', 'execution.time')

#
# SIMLR
context("SIMLR")
normal <- SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0)
if.impute <- SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, if.impute = TRUE)
normalise <- SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0, normalize = TRUE)
test_that("structure of output is compliant", {
  expect_true(all(expected_names %in% names(normal)))
  expect_true(all(expected_names %in% names(if.impute)))
  expect_true(all(expected_names %in% names(normalise)))
})

#
# CIMLR
context("CIMLR")
cimlr <- CIMLR(X = GliomasReduced$in_X, c = 3, cores.ratio = 0)
test_that("structure of output is compliant", {
  expect_true(all(expected_names %in% names(cimlr)))
})

#
# Ranking
context("SIMLR ranking")
ranks <- SIMLR_Feature_Ranking(A = BuettnerFlorian$results$S, X = BuettnerFlorian$in_X)
test_that("structure of output is compliant", {
  expect_equal(names(ranks), c("pval", "aggR"))
})
