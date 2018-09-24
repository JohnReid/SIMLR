#!/usr/bin/env Rscript

#
# Run SIMLR on the neuronal data set from Usoskin et al.
#

#
# load the SIMLR R package
# devtools::document('../..')
# devtools::build('../..')
devtools::load_all('../..')
ls('package:SIMLR')

#
# Configure
set.seed(11111)
.data <- Usoskin
data_set <- 'Usoskin'

#
# Run and summarise SIMLR
res <- run_SIMLR(.data = .data)
summarise_SIMLR(res = res, .data = .data, data_set = data_set)
