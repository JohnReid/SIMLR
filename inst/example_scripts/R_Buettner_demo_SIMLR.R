#!/usr/bin/env Rscript

#
# Run SIMLR on the Buettner data set
#

#
# Load packages
library(tidyverse)

#
# load the SIMLR R package
devtools::load_all('../..')
# devtools::build('../..')
# devtools::document('../..')

#
# run SIMLR
set.seed(11111)
cat("Performing analysis for Buettner data","\n")
.data <- BuettnerFlorian
output_file <- purrr::partial(file.path, 'output')
dir.create(output_file('.', showWarnings = FALSE))
res = SIMLR(X = .data$in_X, c = .data$n_clust, return_intermediaries = TRUE)
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
# make the scatter plots
pdf(output_file('Buettner-scatter.pdf'), width=9, height=6, paper='special')
plot(res$ydata,
     col = get_palette(3)[.data$true_labs[,1]],
     xlab = "SIMLR component 1",
     ylab = "SIMLR component 2",
     pch = 20,
     main = "SIMILR 2D visualization for Buettner data set")
dev.off()


#
# Make a heatmap of S
devtools::load_all('../..')
pdf(output_file('Buettner-heatmap.pdf'), width=8, height=8, paper='special')
similarity.heatmap(res$S,
                   label = str_c('label ', .data$true_labs[,1]),
                   cluster = str_c('cluster ', res$y$cluster))
dev.off()


#
# Show and save timings
print(res$timings)
readr::write_csv(res$timings, output_file('Buettner-timings.csv'))

#
# Show NMI
message('NMI: ', nmi_1)
