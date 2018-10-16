
if (!require("parallel")) install.packages("parallel")
library("parallel")

"FZINB.matrix.fork" <- function(X_counts, Theta0 = NULL, truncate.ratio = 0.3, n_gene = 1000 , cores.ratio = 1){
  if (is.null(Theta0)){
    Theta0 <- prior.zinb(X_counts)$Theta0
    if(Theta0[2]<0){
      warning("A Priori parameter estimate of ZINB: Negative size parameter. Reset to 1.")
      Theta0[2] <- 1
    }
  }
  n_cell <- ncol(X_counts)
  cat(paste("-- Computing MLE for fitting ZINB.\n" ))
  THETA <- MLE.zinb(X_counts,Theta0)
  extract <- which(THETA[,3]>truncate.ratio & THETA[,3]<0.9)
  
  extract_sorted <- head(sort(apply(X_counts[extract,], 1, var), decreasing = TRUE,index.return = TRUE)$ix , n_gene)
  sorted <- extract[extract_sorted]
  cat("-- Computing FZINB for Gene.\n")
  # setup a parallelized estimation of the kernels
  FZINB <-list()
  wd_ = getwd()
  
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterExport(cl=cl, varlist = c("sorted","THETA","X_counts","wd_"), envir=environment())
  FZINB <- parLapply(cl, sorted, function(l){
    setwd(wd_)
    source("../ZINB/functions.r")
    I <- I.zinb(THETA[l,])
    tempK <- as.numeric()
    for (j in 1:ncol(X_counts)){
      tempK <- c(tempK, Fisher.zinb(cbind(X_counts[l,],X_counts[l,j]), THETA[l,], I=I))
    }
    return(tempK) 
  })
  stopCluster(cl)
  FZINB <- Reduce("+", FZINB)/n_gene
  return(matrix(FZINB,ncol = n_cell))
}




"FZINB.matrix.sockets" <- function(X_counts, THETA, Theta0 = NULL, truncate.ratio = 0.3, n_gene = 1000, cores.ratio = 1){
  if (is.null(Theta0)){
    Theta0 <- prior.zinb(X_counts)$Theta0
    if(Theta0[2]<0){
      warning("A Priori parameter estimate of ZINB: Negative size parameter. Reset to 1.")
      Theta0[2] <- 1
    }
  }
  n_cell <- ncol(X_counts)
  cat(paste("-- Computing MLE for fitting ZINB.\n" ))
  #THETA <- MLE.zinb(X_counts,Theta0)
  extract <- which(THETA[,3]>truncate.ratio & THETA[,3]<0.9)
  
  extract_sorted <- head(sort(apply(X_counts[extract,], 1, var), decreasing = TRUE,index.return = TRUE)$ix , n_gene)
  sorted <- extract[extract_sorted]
  tempK <- as.numeric()
  Fisher <- array(0,c(n_cell,n_cell)) #Fisher matrix
  
  cat("-- Computing FZINB parallelly for Genes.\n")
  cores = as.integer(cores.ratio * (detectCores() - 2))
  if (cores < 1) {
    cores = 1
  }
  iterGene <- function(l){
    I <- I.zinb(THETA[l,])
    for (j in 1:n_cell){
      tempK <- c(tempK,Fisher.zinb(cbind(X_counts[l,],X_counts[l,j]), THETA[l,], I=I))
    }
    return(tempK) 
  }
  FZINB <- Reduce("+", mclapply(sorted, FUN = iterGene , mc.cores = cores))/n_gene

  return(matrix(FZINB,ncol = n_cell))
}




"FZINB.matrix" <- function(X_counts, Theta0 = NULL, truncate.ratio = 0.3, n_gene = 1000, cores.ratio = NA){
  if (is.null(Theta0)){
    Theta0 <- prior.zinb(X_counts)$Theta0
    if(Theta0[2]<0){
      warning("A Priori parameter estimate of ZINB: Negative size parameter. Reset to 1.")
      Theta0[2] <- 1
    }
  }
  n_cell <- ncol(X_counts)
  cat(paste("-- Computing MLE for fitting ZINB.\n" ))
  THETA <- MLE.zinb(X_counts,Theta0)
  extract <- which(THETA[,3]>truncate.ratio & THETA[,3]<0.9)
  
  extract_sorted <- head(sort(apply(X_counts[extract,], 1, var), decreasing = TRUE,index.return = TRUE)$ix , n_gene)
  sorted <- extract[extract_sorted]
  tempK <- array(0,c(n_cell,n_cell))
  FZINB <- array(0,c(n_cell,n_cell)) #FZINB matrix
  cat("-- Computing FZINB for Gene:\n")
  for (l in sorted){
    I <- I.zinb(THETA[l,])
    for (j in 1:n_cell){
      tempK[,j] <- Fisher.zinb(cbind(X_counts[l,],X_counts[l,j]), THETA[l,], I=I)
    }
    FZINB <- FZINB + tempK
    if (which(sorted==l)%%100==0){
      cat(paste("-- Gene",which(sorted==l),"\n"))
    }
  }
  return(FZINB/n_gene)
}    




X_counts <- round(10^Test_1_mECS$in_X-1.0)
THETA <- MLE.zinb(X_counts,Theta0)
system.time(FZINB <- FZINB.matrix(X_counts, THETA, cores= 4))



