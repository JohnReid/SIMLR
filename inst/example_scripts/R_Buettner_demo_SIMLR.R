#!/usr/bin/env Rscript

#
# Run SIMLR on the Buettner data set
#

library(tidyverse)

#
# load the SIMLR R package
devtools::load_all('../..')

#
# run SIMLR
set.seed(11111)
cat("Performing analysis for Buettner data","\n")
.data <- BuettnerFlorian
output_file <- purrr::partial(file.path, 'output')
dir.create(output_file('.', showWarnings = FALSE))
res = SIMLR(X = .data$in_X, c = .data$n_clust)
# Calculate NMI
nmi_1 = igraph::compare(.data$true_labs[,1], res$y$cluster, method = "nmi")

#
# Report results
print(res$execution.time)
print('Convergence:')
print(res$converge)
print('Weights:')
print(res$alphaK)
# print('k-means:')
# print(res$y)

#
# make the scatter plots
pdf(output_file('Buettner-scatter.pdf'), width=12, height=8, paper='special')
plot(res$ydata,
     col=c(topo.colors(.data$n_clust))[.data$true_labs[,1]],
     xlab="SIMLR component 1",
     ylab="SIMLR component 2",
     pch=20,
     main="SIMILR 2D visualization for Buettner data set")
dev.off()

#
# Show and save timings
print(res$timings)
readr::write_csv(res$timings, output_file('Buettner-timings.csv'))

#
# Show NMI
message('NMI: ', nmi_1)
