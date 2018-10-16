if (!require("ZIM")) install.packages("ZIM")
library("ZIM")
if (!require("devtools")) install.packages("devtools")
library("devtools")
if (!require("lbfgs")) install.packages("lbfgs")
library("lbfgs")
##### setting #####

#y <- rzinb(10000, theta, mu , pi)

####################Three parameter (mu, theta, pi) ##########################




##### score ###
"score.zinb" <- function(y , Theta){
  mu <- Theta[1]
  theta <- Theta[2] 
  pi <- Theta[3]
  U  <- array(0,c(3,length(y)))
  frac <- theta/(theta+mu)
  
  U[1,] <- y/mu - (y+theta)/(mu+theta)
  U[2,] <- digamma(y+theta) - digamma(theta) + log(theta/(theta+mu)) + (mu-y)/(theta+mu)
  U[3,] <- -1/(1-pi)
  
  f0 <- pi + (1-pi)*frac^theta
  U[1,y==0] <- -(1-pi)*frac^(theta+1)/f0
  U[2,y==0] <- (1-pi)*frac^theta*(log(frac)+1-frac)/f0
  U[3,y==0] <- (1-frac^theta)/f0
  
  return(U)
} 

#### Fisher information matrix ####




"I.zinb" <- function(Theta){
  mu <- Theta[1]
  theta <- Theta[2]
  pi <- Theta[3]
  
  I <- matrix(0,nrow=3,ncol=3)
  frac <- 1/(theta+mu)
  
  fNB0 <- dnbinom(0,mu = mu, size = theta)
  fZINB0 <- dzinb(0, lambda = mu, k = theta, omega = pi)
  score.NB0 <- c(-theta*frac , log(theta*frac)+mu*frac)
  
  Hessian.logNB0 <- matrix(c(theta*frac^2,-mu*frac^2,-mu*frac^2,mu*frac/theta-mu*frac^2), nrow=2)
  
  Nabla.NB0  <- fNB0 * score.NB0
  Hessian.NB0 <-  fNB0* (score.NB0%*%t(score.NB0) + Hessian.logNB0)
  
  PhiPhi.logZINB0 <- (1-pi)*(Hessian.NB0*fZINB0 - (1-pi)*Nabla.NB0%*%t(Nabla.NB0))/fZINB0^2
  Phipi.logZINB0 <- (-Nabla.NB0*fZINB0 - (1-pi)*Nabla.NB0*(1-fNB0))/fZINB0^2
  
  n=0
  infisum =0
  add <- (1-pnbinom(n, mu=mu,size = theta))/(theta+n)^2
  while (add > 1e-9){
    add <- (1-pnbinom(n, mu=mu,size = theta))/(theta+n)^2
    infisum <- infisum + add
    n=n+1
  }
  
  I_NB <- diag(c(1/mu -frac, -mu*frac/theta+infisum))
  
  
  I[3,3] <- (1 - fNB0)^2/fZINB0 + (1-fZINB0)/(1-pi)^2
  I[1:2,1:2] <- -PhiPhi.logZINB0*fZINB0 + (1-pi)*(I_NB + Hessian.logNB0*fNB0)
  I[1:2,3] <- -Phipi.logZINB0*fZINB0
  I[3,1:2] <- I[1:2,3]
  
  return(I)
}


"Fisher.zinb" <- function(Y,Theta,iter = 10000, I = NULL){
  if (is.null(I) ) {I <- I.zinb(Theta)}
  L <- t(chol(I))
  if (is.null(dim(Y))){
    Fisher <- t(solve(L,score.zinb(Y[1],Theta)))%*%solve(L,score.zinb(Y[2],Theta))
  }
  else{
    score <- array(0,c(3,dim(Y)[1],2))
    score[,,1] <- score.zinb(Y[,1],Theta)
    score[,,2] <- score.zinb(Y[,2],Theta)
    Fisher <- colSums(solve(L,score[,,1])*solve(L,score[,,2]))
  }
  return(Fisher)
}


####################Two parameter (mu, theta) ##########################
##### score ###
"score.zinb2" <- function(y , Theta){
  mu <- Theta[1]
  theta <- Theta[2] 
  pi <- Theta[3]
  U  <- array(0,c(2,length(y)))
  frac <- theta/(theta+mu)
  
  U[1,] <- y/mu - (y+theta)/(mu+theta)
  U[2,] <- digamma(y+theta) - digamma(theta) + log(theta/(theta+mu)) + (mu-y)/(theta+mu)
  
  
  f0 <- pi + (1-pi)*frac^theta
  U[1,y==0] <- -(1-pi)*frac^(theta+1)/f0
  U[2,y==0] <- (1-pi)*frac^theta*(log(frac)+1-frac)/f0
  
  
  return(U)
} 

"I.zinb2" <- function(Theta){
  return(I.zinb(Theta)[1:2,1:2])
}

"Fisher.zinb2" <- function(Y,Theta,iter = 10000, I = NULL){
  mu <- Theta[1]
  theta <- Theta[2] 
  pi <- Theta[3]
  if (is.null(I) ) {I <- I.zinb2(Theta)}
  L <- t(chol(I))
  if (is.null(dim(Y))){
    Fisher <- t(solve(L,score.zinb2(Y[1],Theta)))%*%solve(L,score.zinb2(Y[2],Theta))
  }
  else{
    score <- array(0,c(2,dim(Y)[1],2))
    score[,,1] <- score.zinb2(Y[,1],Theta)
    score[,,2] <- score.zinb2(Y[,2],Theta)
    Fisher <- colSums(solve(L,score[,,1])*solve(L,score[,,2]))
  }
  return(Fisher)
}

#### parameter estimation in ZINB #####

"prior.zinb" <- function(X_counts){
  
  mu0 <- apply(X_counts,1, FUN = function(x) mean(x[x!=0]))
  theta0 <- mu0^2/(apply(X_counts, 1, var)-mu0)
  pi0 <- rowSums(X_counts==0)/dim(X_counts)[2]
  Theta0_mean <- c(mean(mu0),mean(theta0),mean(pi0))
  Theta0_var <- c(var(mu0),var(theta0),var(pi0))
  
  result <- list()
  result[["Theta0"]] = Theta0_mean
  result[["variance"]] = Theta0_var
  
  return(result)
}


###### bfgs in parameter mu, theta, pi #######

## INPUT : matrix with genes x cells
## OUTPUT: matrix with genes x 3, each row represents the MLE estimate for that gene


"MLE.zinb" <- function(X,Theta, invisible = 1){
  output <- numeric()
  if (is.null(dim(X))){
    X <- t(matrix(X))
  }
  mu0 <- Theta[1]
  theta0 <- Theta[2] #aggregation paramter, dispersion parameter
  pi0 <- Theta[3]
  phi <- c(log(mu0),log(theta0),qnorm(pi0))
  for (l in 1:dim(X)[1]){
    
    dat <- X[l,]
    
    objective <- function(phi){
      mu <- exp(phi[1])
      theta <- exp(phi[2])
      pi <- pnorm(phi[3])
      n0 <- sum(dat==0)
      datt <- round(dat[dat>0])
      n <- length(dat)
      
      f<- n0 * log(pi + (1-pi)*(theta/(theta+mu))^theta) +
        (n-n0)*log(1-pi) + sum(log(dnbinom(datt,size = theta, mu = mu )))
      return(-f)   
    }
    
    
    gradient <- function(phi){
      mu <- exp(phi[1])
      theta <- exp(phi[2])
      pi <- pnorm(phi[3])
      n0 <- sum(dat==0)  
      datt <- dat[dat>0]
      n <- length(dat)
      
      ratio <- theta/(theta+mu)
      dmu <- -n0*(1-pi)*
        ratio^(theta+1)/(pi + (1-pi)*ratio^theta)+
        sum(datt/mu - (datt+theta)/(mu+theta))
      
      
      dtheta <- n0* (1-pi) *ratio^theta *(log(ratio)+1-ratio)/(pi+(1-pi)*ratio^theta) +  
        sum(digamma(datt+theta) - digamma(theta) + log(ratio) + (mu-datt)/(theta+mu))
      
      dpi <- n0*(1-ratio^theta)/(pi + (1-pi)*ratio^theta)-
        (n-n0)/(1-pi)
      
      dphi1 <- dmu*exp(phi[1])
      dphi2 <- dtheta*exp(phi[2])
      dphi3 <- dpi*dnorm(phi[3])
      
      return(c(-dphi1,-dphi2,-dphi3))
    }
    
    out<- lbfgs(objective, gradient, phi, invisible = invisible)
    #out <- optim(par = phi, fn = objective, gr = gradient ,method = "BFGS" )   
    output <- rbind(output,c(exp(out$par[1]),
                             exp(out$par[2]),
                             pnorm(out$par[3])))
  }
  colnames(output) <- c("mu","theta","pi")
  return(output)
}


###FZINB matrix ####
# more other version of FZINB.matrix function, such as non-parallelling version, can be seen in file tryParallel.R
"FZINB.matrix" <- function(X_counts, Theta0 = NULL, truncate.ratio = 0.3, n_gene = 1000 , cores.ratio = 1){
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
    source("./functions.r")
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


"FZINB.Distance.matrix" <- function(FZINB){
  k <- 1/sqrt(diag(FZINB))
  FD_matrix <- sqrt(abs(2 - 2*FZINB * (k %*% t(k))))
  diag(FD_matrix)<-0
  return(FD_matrix)
} 
