#
# A file to store miscellaneous debugging code
#

library(microbenchmark)
library(ggplot2)
devtools::load_all('../..')


devtools::test(filter = 'alt_calcs')

#
# Retrieve one of the S intermediaries
S <- res$intermediaries$S[5, , ]
dim(S)
all.equal(rowSums(S), colSums(S))
summary(rowSums(S))
S[1:5, 1:5]
mean(S != 0)

#
# network.diffusion parameters
A <- S
K <- 10
scale_by_DD <- TRUE


devtools::load_all('../..')
L.orig <- laplacian(S, 'orig')
L.apply <- laplacian(S, 'apply')
L.orig[1:5, 1:5]
L.apply[1:5, 1:5]
all.equal(laplacian(S, 'orig'), laplacian(S, 'apply'))
microbenchmark(laplacian(S, 'orig'), laplacian(S, 'apply'))


#
# Look at eigenvalue adjustment
df. <- data.frame(original = d, transformed = ((1 - alpha) * d) / (1 - alpha * d^beta), i = 1:length(d))
df.m <- reshape2::melt(df., id.vars = 'i', variable.name = 'eigenvalue', value.name = 'd')
gp.by.dim <- ggplot(df.m, aes(x = i, y = d, colour = eigenvalue)) + geom_point()
gp.scatter <- ggplot(df., aes(x = original, y = transformed)) + geom_point()
gridExtra::grid.arrange(gp.by.dim, gp.scatter, nrow = 2)

#
#
microbenchmark(divide_rows_by_diag(W, method = 'sweep'),
               divide_rows_by_diag(W, method = 'apply'))


#
# Check how to write diag(v) %*% X as a sweep
nr <- 511
nc <- 513
X <- matrix(rnorm(nr * nc), nrow = nr)
v <- rnorm(nr)
fast <- sweep(X, MARGIN = 1, STATS = v, FUN = '*')
all.equal(multiply_rows(v, X, 'sweep'), multiply_rows(v, X, 'diag'))
microbenchmark(multiply_rows(v, X, 'sweep'), multiply_rows(v, X, 'diag'))
v <- rnorm(nc)
all.equal(multiply_cols(v, X, 'sweep'), multiply_cols(v, X, 'diag'))
microbenchmark(multiply_cols(v, X, 'sweep'), multiply_cols(v, X, 'diag'))


