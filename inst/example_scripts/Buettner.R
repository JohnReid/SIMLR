#!/usr/bin/env Rscript

#
# Run SIMLR on the mouse embryonic stem cell data set from Buettner et al.
#

library(tidyverse)
library(ggthemes)

theme_set(theme_few())
scale_colour_discrete <- function(...) scale_colour_few()
scale_fill_discrete <- function(...) scale_fill_few()

#
# load the SIMLR R package
# devtools::document('../..')
# devtools::build('../..')
devtools::load_all('../..')
# ls('package:SIMLR')

#
# Configure
set.seed(11111)
.data <- Buettner
data_set <- 'Buettner'

#
# Estimate number of clusters
# num_clust_estimate <- SIMLR_Estimate_Number_of_Clusters(.data$in_X)

#
# Run and summarise SIMLR
res <- SIMLR::run_SIMLR(.data = .data)
SIMLR::summarise_SIMLR(res = res, .data = .data, data_set = data_set)



A <- res$S
K <- 10
diag(A) = 0
P = dominate.set(abs(A), min(K, nrow(A) - 1)) * sign(A)
mean(P == 0)  # 93% sparse for small example
A[1:6,1:6]
P[1:6,1:6]

# sum the absolute value of each row of P
DD = apply(abs(P), MARGIN = 1, FUN = sum)

# set the diagonal of P to be DD + 1
diag(P) = DD + 1

# compute the transition field of P
P = transition.fields(P)

# compute the eigenvalues and eigenvectors of P
eigen_P = eigen(P)
U = eigen_P$vectors
D = eigen_P$values

# set to d the real part of the diagonal of D
d = Re(D + .Machine$double.eps)

# perform the diffusion
alpha = 0.8
beta = 2
d = ((1-alpha)*d)/(1-alpha*d^beta)

# set D to be a diagonal matrix of the real part of d
D = diag(Re(d))

# finally compute W
W = U %*% D %*% t(U)
diagonal_matrix = diag(rep(1, nrow(W)))
W = (W * (1 - diagonal_matrix)) / apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=(1-diag(W))})
if( scale_by_DD ) {
  # This line is missing in network.diffusion.numc()
  W = diag(DD) %*% W
}
# Ensure W symmetric
W = (W + t(W)) / 2
# Ensure all W are non-negative
W[W < 0] = 0

