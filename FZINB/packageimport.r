if (!require("Matrix")) install.packages("Matrix")
library(Matrix)
if (!require("parallel")) install.packages("parallel")
library(parallel)

# load the igraph package to compute the NMI
if (!require("igraph")) install.packages("igraph")
library(igraph)

# load the palettes for the plots
if (!require("grDevices")) install.packages("grDevices")
library(grDevices)



# load the SIMLR R package
source("../R/SIMLR.R")
source("../R/compute.multiple.kernel.R")
source("../R/network.diffusion.R")
source("../R/utils.simlr.R")
source("../R/tsne.R")
source("../R/calc.DD.R")

source("./functions.r")

source("./SIMLR_FZINB.r")
source("./compute.Gaussian.kernel.r")

