#!/usr/bin/env Rscript

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
.data <- ZeiselAmit
data_set <- 'ZeiselAmit-LS'

#
# Estimate number of clusters
# num_clust_estimate <- SIMLR_Estimate_Number_of_Clusters(.data$in_X)

#
# Subsample cells
# cells <- sample(ncol(.data$in_X), 340)
# .data$in_X <- .data$in_X[, cells]
# .data$true_labs <- as.matrix(.data$true_labs[cells, ], ncol=1)

#
# Run and summarise SIMLR (large scale)
resLS <- SIMLR::SIMLR_Large_Scale(X = .data$in_X, c = .data$n_clust)
SIMLR::summarise_SIMLR(res = resLS, .data = .data, data_set = data_set)

#
# Run and summarise SIMLR
dim(res$S)
dim(res$intermediaries$S)
nrow(res$intermediaries$S[13, 1:5, 1:5])
res <- SIMLR::run_SIMLR(.data = .data)

devtools::load_all('../..')
SIMLR::summarise_SIMLR(res = res, .data = .data, data_set = 'ZeiselAmit')

traceback()
