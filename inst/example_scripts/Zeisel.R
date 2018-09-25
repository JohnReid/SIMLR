#!/usr/bin/env Rscript

#
# Run SIMLR on the mouse brain data set from Zeisel et al.
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
.data <- Zeisel
data_set <- 'Zeisel'

#
# Run and summarise SIMLR
res <- SIMLR::run_SIMLR(.data = .data)
SIMLR::summarise_SIMLR(res = res, .data = .data, data_set = data_set)
