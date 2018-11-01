context("CIMLR")
library(SIMLR)

test_that("CIMLR runs without errors", {
  set.seed(11111)
  CIMLR(X = GliomasReduced$in_X, c = 3, cores.ratio = 0)
})
