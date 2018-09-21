#!/usr/bin/env Rscript

#
# Run SIMLR on the Buettner data set
#


#
# load the SIMLR R package
devtools::load_all('../..')
# devtools::build('../..')
# devtools::document('../..')

#
# test SIMLR.R on example 1
set.seed(11111)
cat("Performing analysis for Buettner data","\n")
.data <- BuettnerFlorian
res = SIMLR(X = .data$in_X, c = .data$n_clust)
print(res$execution.time)
print('Convergence:')
print(res$converge)
print('Weights:')
print(res$alphaK)
# print('k-means:')
# print(res$y)
nmi_1 = igraph::compare(.data$true_labs[,1], res$y$cluster, method="nmi")

#
# make the scatter plots
dir.create('tests', showWarnings = FALSE)
pdf('tests/Buettner-scatter.pdf', width=12, height=8, paper='special')
plot(res$ydata,
     col=c(topo.colors(.data$n_clust))[.data$true_labs[,1]],
     xlab="SIMLR component 1",
     ylab="SIMLR component 2",
     pch=20,
     main="SIMILR 2D visualization for Buettner data set")
dev.off()

#
# Show timings
res$timings

#
# Show NMI
message('NMI: ', nmi_1)
