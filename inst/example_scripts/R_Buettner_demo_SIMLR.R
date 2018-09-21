#!/usr/bin/env Rscript

#
# Run SIMLR on the Buettner data set
#

# required external packages for SIMLR
# library(Matrix)
# library(parallel)

# load the igraph package to compute the NMI
# library(igraph)

# load the palettes for the plots
# library(grDevices)

# load the SIMLR R package
devtools::load_all('../..')



# load the C file

# NOTE 1: we make use of an external C program during the computations of SIMLR.
# The code is located in the R directory in the file projsplx_R.c. In order to 
# use SIMLR one needs to compite the program. To do so, one needs to run on the 
# shell the command R CMD SHLIB -c projsplx_R.c. 
# The precompiled projsplx_R.so is already provided for MAC OS X only. 
# If one wants to use SIMLR on other operative systems, the file projsplx_R.so 
# needs to be deleted, and re-compiled. 

# NOTE 2: for Windows, the command dyn.load("./R/projsplx_R.so") needs to be 
# substituted with the command dyn.load("./R/projsplx_R.dll"). 

dyn.load("./R/projsplx_R.so")

# load the datasets
load(file="./data/Test_1_mECS.RData")

# test SIMLR.R on example 1
set.seed(11111)
cat("Performing analysis for Test_1_mECS","\n")
res_example1 = SIMLR(X=Test_1_mECS$in_X,c=Test_1_mECS$n_clust)
print(res_example1$execution.time)
print('Convergence:')
print(res_example1$converge)
print('Weights:')
print(res_example1$alphaK)
# print('k-means:')
# print(res_example1$y)
nmi_1 = compare(Test_1_mECS$true_labs[,1],res_example1$y$cluster,method="nmi")

# make the scatter plots
dir.create('tests', showWarnings = FALSE)
pdf('tests/Buettner-scatter.pdf', width=12, height=8, paper='special')
plot(res_example1$ydata,
     col=c(topo.colors(Test_1_mECS$n_clust))[Test_1_mECS$true_labs[,1]],
     xlab="SIMLR component 1",
     ylab="SIMLR component 2",
     pch=20,
     main="SIMILR 2D visualization for Test_1_mECS")
dev.off()

#
# Show timings
res_example1$timings

#
# Show NMI
message('NMI: ', nmi_1)
