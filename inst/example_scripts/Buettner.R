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
# Run SIMLR
set.seed(11111)
res <- apply_SIMLR(Buettner, 'Buettner')
