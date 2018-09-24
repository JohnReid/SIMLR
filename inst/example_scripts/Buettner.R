#!/usr/bin/env Rscript

#
# Run SIMLR on the mouse embryonic stem cell data set from Buettner et al.
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
.data <- Buettner
data_set <- 'Buettner'

#
# Run and summarise SIMLR
res <- SIMLR::run_SIMLR(.data = .data)
SIMLR::summarise_SIMLR(res = res, .data = .data, data_set = data_set)
