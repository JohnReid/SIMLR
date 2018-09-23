#!/usr/bin/env Rscript

#
# Run SIMLR on the neuronal data set from Usoskin et al.
#

#
# Load packages
library(tidyverse)

#
# load the SIMLR R package
# devtools::document('../..')
# devtools::build('../..')
devtools::load_all('../..')
ls('package:SIMLR')

#
# Configuration
.data <- Usoskin
output_file <- purrr::partial(file.path, 'output', 'Usoskin')
dir.create(output_file('.'), recursive = TRUE, showWarnings = FALSE)

#
# run SIMLR
set.seed(11111)
message("Running SIMLR")
res = SIMLR(X = .data$in_X, c = .data$n_clust, return_intermediaries = TRUE)

#
# Calculate NMI
nmi_1 = igraph::compare(.data$true_labs[,1], res$y$cluster, method = "nmi")

#
# Report results
print('Iterations:')
print(res$iter)
print(res$execution.time)
print('Convergence:')
print(res$converge)
print('Weights:')
print(res$alphaK)

#
# Scatter plot of dimensionality reduction
pdf(output_file('scatter.pdf'), width=9, height=6, paper='special')
plot(res$ydata,
     col = get_palette(3)[.data$true_labs[,1]],
     xlab = "SIMLR component 1",
     ylab = "SIMLR component 2",
     pch = 20,
     main = "SIMILR 2D visualization")
dev.off()


#
# Make a heatmap of S
# devtools::load_all('../..')
pdf(output_file('heatmap.pdf'), width=8, height=8, paper='special')
similarity.heatmap(res$S,
                   label = str_c('label ', .data$true_labs[,1]),
                   cluster = str_c('cluster ', res$y$cluster))
dev.off()

#
# Show and save timings
print(res$timings)
readr::write_csv(res$timings, output_file('timings.csv'))

#
# Show NMI
message('NMI: ', nmi_1)
