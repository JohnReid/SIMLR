#!/usr/bin/env Rscript

#
# Run SIMLR on the cerebral cortex data set from Pollen et al.
#

#
# load the SIMLR R package
# devtools::document('../..')
# devtools::build('../..')
devtools::load_all('../..')
# ls('package:SIMLR')

#
# Configure
set.seed(11111)
.data <- Pollen
data_set <- 'Pollen'

#
# Run and summarise SIMLR
res <- SIMLR::run_SIMLR(.data = .data)
SIMLR::summarise_SIMLR(res = res, .data = .data, data_set = data_set)
