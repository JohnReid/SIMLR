#!/usr/bin/env Rscript

#
# A script to examine the internals of SIMLR for the purpose of understanding its implementation
#

# required external packages for SIMLR
library(Matrix)
library(parallel)

# load the igraph package to compute the NMI
library(igraph)

# load the palettes for the plots
library(grDevices)

# For plotting
library(ggplot2)

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
# Plot distances
pdf('tests/distances.pdf')
qplot(as.vector(Diff[upper.tri(Diff)]))
dev.off()
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
D_Kernels = multiple.unnorm.kernels(t(X), cores.ratio=1)
idxs = sample(length(D_Kernels), size=3)
for (i in idxs) {
  print(D_Kernels[[i]][1:4, 1:4])
}
D_Dists = norm.and.calc.dists(D_Kernels)
for (i in idxs) {
  print(D_Dists[[i]][1:4, 1:4])
}
# saveRDS(D_Dists, 'tests/D_kernels.rds')
class(D_Dists)
names(D_Dists)
length(D_Dists)
dim(D_Dists[[1]])
G = D_Kernels[[idxs[1]]]
dists = D_Dists[[idxs[1]]]
uptri = upper.tri(G)
.df = data.frame(
  k = as.vector(G[uptri]),
  d = as.vector(dists[uptri]))
pdf('tests/kernel-vs-dist.pdf')
ggplot(.df, aes(x = k, y = d)) + geom_point(alpha = .1)
dev.off()


#
# set up some parameters
#
NITER = 30
num = ncol(X)
r = -1
beta = 0.8
k = 10
c = 3  # number of clusters
# alphaK looks like a Dirichlet prior: 1 / # categories
alphaK = 1 / rep(length(D_Kernels), length(D_Kernels))
#
# distX is the average of the distances
distX = Reduce("+", D_Kernels) / length(D_Kernels)
#
# sort each row of distX into distX1 and retain the ordering vectors in idx
res = apply(distX, MARGIN = 1, FUN = function(x) return(sort(x, index.return = TRUE)))
distX1 = array(0, c(nrow(distX), ncol(distX)))
idx = array(0, c(nrow(distX), ncol(distX)))
for(i in 1:nrow(distX)) {
    distX1[i,] = res[[i]]$x
    idx[i,] = res[[i]]$ix
}

# Unsure what A is so far
A = array(0, c(num, num))
# di contains the distances to the (k+1) nearest neighbours
di = distX1[,2:(k+2)]
# rr is half the difference between k times the k+1'th nearest neighbour distance and
# the sum of the k nearest neighbour distances
rr = 0.5 * (k * di[, k + 1] - apply(di[, 1:k], MARGIN = 1, FUN = sum))
# Index the k+1 nearest neighbours
id = idx[, 2:(k + 2)]
# The numerator 
dim(di)
k
numerator = apply(array(0, c(length(di[,k+1]), dim(di)[2])),
                  MARGIN = 2,
                  FUN = function(x) {x=di[, k+1]}) - di

numerator.new = t(di[, k+1] -t(di))
all.equal(numerator, numerator.new)

temp = (k*di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum) + .Machine$double.eps)
denominator = apply(array(0,c(length(temp),dim(di)[2])),MARGIN=2,FUN=function(x) {x=temp})
temp = numerator / denominator
a = apply(array(0,c(length(t(1:num)),dim(di)[2])),MARGIN=2,FUN=function(x) {x=1:num})
A[cbind(as.vector(a),as.vector(id))] = as.vector(temp)
if(r<=0) {
    r = mean(rr)
}
lambda = max(mean(rr),0)
A[is.nan(A)] = 0
A0 = (A + t(A)) / 2
S0 = max(max(distX)) - distX

cat("Performing network diffiusion.\n")

# perform network diffiusion
S0 = network.diffusion(S0,k)

#
# Test dn
w = S0
D = apply(w, MARGIN=2, FUN=sum)
#
# type "ave" returns D^-1*W
D = 1 / D
D_temp = matrix(0, nrow=length(D), ncol=length(D))
D_temp[cbind(1:length(D),1:length(D))] = D
D_temp.new = diag(D)
all.equal(D_temp, D_temp.new)
wn = D_temp %*% w
rowSums(wn)
colSums(wn)

# Normalise S0
S0.gph = dn(S0, 'gph')
S0 = dn(S0, 'ave')
all.equal(S0, S0.gph)
all.equal(t(S0.gph), S0.gph)
S = S0
D0 = diag(apply(S,MARGIN=2,FUN=sum))
L0 = D0 - S

eig1_res = eig1(L0, c, 0)
F_eig1 = eig1_res$eigvec
temp_eig1 = eig1_res$eigval
evs_eig1 = eig1_res$eigval_full
