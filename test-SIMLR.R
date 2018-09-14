# required external packages for SIMLR
library(Matrix)
library(parallel)

# load the igraph package to compute the NMI
library(igraph)

# load the palettes for the plots
library(grDevices)

# load the SIMLR R package
source("./R/SIMLR.R")
source("./R/compute.multiple.kernel.R")
source("./R/network.diffusion.R")
source("./R/utils.simlr.R")
source("./R/tsne.R")

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
# load(file="./data/Test_2_Kolod.RData")
# load(file="./data/Test_3_Pollen.RData")
# load(file="./data/Test_4_Usoskin.RData")

#
# Check data
dim(Test_1_mECS$in_X)  # genes x samples
Test_1_mECS$n_clust  # number of clusters
X <- Test_1_mECS$in_X
x <- t(X)
#
# Parameters
sigma = seq(2,1,-0.25)
allk = seq(10,30,2)
#
# Calculate distances
Diff = dist2(x)^2  # Distances squared square (i.e. to power 4) between samples
Diff[1:3,1:3]
#
# Check distances
dim(x)
# Calculate the distances squared
d12sq <- sum((x[1,] - x[2,])^2)
d13sq <- sum((x[1,] - x[3,])^2)
d23sq <- sum((x[2,] - x[3,])^2)
# Check the distances squared match the Diff
d12sq^2 - Diff[1, 2]
d13sq^2 - Diff[1, 3]
d23sq^2 - Diff[2, 3]
#
# Check sorting
Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))
Diff_sort[1:3,1:3]
Diff_sort_fun <- Diff_sort
#
# Mean of k-nearest-neighbours
k <- 10
TT = apply(Diff_sort_fun[,2:(k+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
# Do an outer average
Sig2 = outer(TT, TT, FUN = function(x, y) (x + y) / 2)
TT = matrix(data = TT, nrow = length(TT), ncol = 1)
dim(TT)
Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
dim(Sig)
Sig[1:3, 1:3]
matrix(rep(TT, 3), ncol=3)[1:3,]
Sig = Sig + t(Sig)
Sig = Sig / 2
max(abs(Sig - Sig2))
all.equal(Sig, Sig2)

#
# Ensure entries are valid, i.e. at least as great as the machine precision!?!?
Sig2 = Sig
Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
Sig = Sig * Sig_valid + .Machine$double.eps
all.equal(Sig, Sig2)

# N.B. Diff_fun == Diff == the squared squared distance (i.e. power of 4)
sigma = 1.5
Sig
W = dnorm(Diff, mean = 0, sd = 1.5*Sig)
K = (W + t(W)) / 2

#
# Test kernel normalisation
k = 1/sqrt(diag(K)+1)  # In practice this is a vector of ones
G = K * (k %*% t(k))   # So k %*% t(k) is a matrix of ones
G1.old = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
G1 = matrix(rep(diag(G), nrow(G)), nrow=nrow(G))
G1[1:4, 1:4]
all.equal(G1, G1.old)
D_Kernels_tmp = (G1 + t(G1) - 2*G)/2
diag(D_Kernels_tmp)
D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
test.idxs <- seq.int(1, 181, length.out=4)
test.idxs <- 1:4
K[test.idxs, test.idxs]
D_Kernels_tmp[test.idxs, test.idxs]
cor(as.vector(K), as.vector(D_Kernels_tmp))  # Negative correlation double check that D is a distance
all.equal(D_Kernels_tmp, kernel.distance.2(kernel.normalise(K)))

#
# Kernels
D_Kernels = multiple.kernel(t(X), cores.ratio=1)
for (i in 1:length(D_Kernels)) {
  print(D_Kernels[[i]][1:4, 1:4])
}
# saveRDS(D_Kernels, 'D_kernels.rds')
class(D_Kernels)
names(D_Kernels)
length(D_Kernels)
dim(D_Kernels[[1]])



# # test SIMLR.R on example 1
set.seed(11111)
cat("Performing analysis for Test_1_mECS","\n")
res_example1 = SIMLR(X=Test_1_mECS$in_X,c=Test_1_mECS$n_clust)
nmi_1 = compare(Test_1_mECS$true_labs[,1],res_example1$y$cluster,method="nmi")

# # test SIMLR.R on example 2
# set.seed(22222)
# cat("Performing analysis for Test_2_Kolod","\n")
# res_example2 = SIMLR(X=Test_2_Kolod$in_X,c=Test_2_Kolod$n_clust)
# nmi_2 = compare(Test_2_Kolod$true_labs[,1],res_example2$y$cluster,method="nmi")

# # test SIMLR.R on example 3
# set.seed(33333)
# cat("Performing analysis for Test_3_Pollen","\n")
# res_example3 = SIMLR(X=Test_3_Pollen$in_X,c=Test_3_Pollen$n_clust)
# nmi_3 = compare(Test_3_Pollen$true_labs[,1],res_example3$y$cluster,method="nmi")

# # test SIMLR.R on example 4
# set.seed(44444)
# cat("Performing analysis for Test_4_Usoskin","\n")
# res_example4 = SIMLR(X=Test_4_Usoskin$in_X,c=Test_4_Usoskin$n_clust)
# nmi_4 = compare(Test_4_Usoskin$true_labs[,1],res_example4$y$cluster,method="nmi")

# make the scatterd plots
pdf('scatter.pdf', width=12, height=8, paper='special')
# par(mfrow=c(2,2))
plot(res_example1$ydata,
     col=c(topo.colors(Test_1_mECS$n_clust))[Test_1_mECS$true_labs[,1]],
     xlab="SIMLR component 1",
     ylab="SIMLR component 2",
     pch=20,
     main="SIMLR 2D visualization for Test_1_mECS")
# plot(res_example2$ydata,
#      col=c(topo.colors(Test_2_Kolod$n_clust))[Test_2_Kolod$true_labs[,1]],
#      xlab="SIMLR component 1",
#      ylab="SIMLR component 2",
#      pch=20,
#      main="SIMLR 2D visualization for Test_2_Kolod")
# plot(res_example3$ydata,
#      col=c(topo.colors(Test_3_Pollen$n_clust))[Test_3_Pollen$true_labs[,1]],
#      xlab="SIMLR component 1",
#      ylab="SIMLR component 2",
#      pch=20,
#      main="SIMLR 2D visualization for Test_3_Pollen")
# plot(res_example4$ydata,
#      col=c(topo.colors(Test_4_Usoskin$n_clust))[Test_4_Usoskin$true_labs[,1]],
#      xlab="SIMLR component 1",
#      ylab="SIMLR component 2",
#      pch=20,
#      main="SIMLR 2D visualization for Test_4_Usoskin")
dev.off()
