
"SIMLR.COMBINE" <- function(X,
                  c,
                  no.dim = NA,
                  k = 10,
                  if.impute = FALSE,
                  normalize = FALSE,
                  cores.ratio = 1,
                  return_intermediaries = FALSE )
{
  # convert SCESet
  if (is(X, "SCESet")) {
    message("X is an SCESet, converting to input matrix.")
    X = X@assayData$exprs
  }

  
  # Measure execution times of sub-tasks
  timer = CumulativeTimer$new()
  
  # set any required parameter to the defaults
  if(is.na(no.dim)) {
    no.dim = c
  }
  
  # check the if.impute parameter
  if(if.impute == TRUE) {
    X = t(X)
    X_zeros = which(X==0,arr.ind=TRUE)
    if(length(X_zeros)>0) {
      R_zeros = as.vector(X_zeros[,"row"])
      C_zeros = as.vector(X_zeros[,"col"])
      ind = (C_zeros - 1) * nrow(X) + R_zeros
      X[ind] = as.vector(colMeans(X))[C_zeros]
    }
    X = t(X)
  }
  timer$add('impute')
  
  # check the normalize parameter
  if(normalize == TRUE) {
    X = t(X)
    X = X - min(as.vector(X))
    X = X / max(as.vector(X))
    C_mean = as.vector(colMeans(X))
    X = apply(X,MARGIN=1,FUN=function(x) return(x-C_mean))
  }
  timer$add('normalise')
  
  
  
  ## involvement of FZINB kernel
  execution.time = list()
  ptm = proc.time()
  cat("Computing the FZINB Kernels.\n")
  X_counts <- round(10^X-1.0)
  FZINB <- FZINB.matrix(X_counts, Theta0 = NULL, n_gene = 1000, cores.ratio = cores.ratio, LogCOUNTS = FALSE, distance = FALSE)
  execution.time["FZINB_kernel"] = (proc.time() - ptm)[1]
  
  ptm = proc.time()
  cat("Computing the Gaussian Kernels based on Euclidean distance space.\n")
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = start_cluster( cores )
  D_Kernels = multiple.kernel(t(X), cl = cl)
  execution.time["Gaussian_kernel"] =(proc.time() - ptm)[1]
  
  cat("Computing the Gaussian Kernels based on FZINB distance space.\n")
  FD_Kernels <- Gaussian.FZINB.kernel(t(X), cores.ratio = cores.ratio, FZINB = FZINB)
  
  cat("Computing the FZINB kernel distance.\n")
  D_FZINB <- Matrix(FZINB.Distance.matrix(FZINB))
  
  results_all <- list()
  methodname <- c("SIMLR","G_FZINBD","FZINB", "Threekind")
  for (Kernels in list("SIMLR" = D_Kernels,
                       "G_FZINBD" = FD_Kernels,
                       "FZINB" = append(D_FZINB,D_FZINB),
                       "Threekind" = append(append(D_Kernels,FD_Kernels),D_FZINB)  )){
    rm(list=setdiff(ls(), c("Kernels",
                            "results_all",
                            "D_Kernels",
                            "FD_Kernels",
                            "D_FZINB",
                            "FZINB",
                            "methodname",
                            "X",
                            "c", 
                            "no.dim",
                            "k",
                            "if.impute",
                            "normalize",
                            "cores.ratio",
                            "execution.time",
                            "return_intermediaries",
                            "timer",
                            "cl") ))

    # set some parameters
    NITER = 30
    num = ncol(X)
    beta = 0.8
    
    # Create arrays to store intermediaries if requested to
    if( return_intermediaries ) {
      intermediaries = list(S = array(NA, dim = c(NITER+1, num, num)),
                            Snd = array(NA, dim = c(NITER+1, num, num)),
                            Lvec =  array(NA, dim = c(NITER+1, c, num)),
                            Lval =  array(NA, dim = c(NITER+1, c)),
                            alphaK = array(NA, dim = c(NITER+1, length(D_Kernels))),
                            dists = array(NA, dim = c(NITER+1, num, num)))
    }
    #
    # alphaK looks like a Dirichlet prior: 1 / # categories
    # that is even weights across the kernels
    alphaK = 1 / rep(length(D_Kernels), length(D_Kernels))
    if( return_intermediaries ) intermediaries$alphaK[1,] = alphaK
    #
    # distX is the average of the distances
    distX = Reduce("+", D_Kernels) / length(D_Kernels)
    if( return_intermediaries ) intermediaries$dists[1,,] = as.matrix(distX)
    timer$add('calc.distances')
    #
    # sort each row of distX into distX1 and retain the ordering vectors in idx
    temp = sort.rows(distX, cl = cl)
    distX1 = temp$sorted
    idx = temp$idx
    timer$add('sort.distances')
    
    #
    # Use data to determine lambda
    #
    # di contains the distances to the (k+1) nearest neighbours
    di = distX1[,2:(k+2)]
    #
    # We don't use distX1 further down
    rm(distX1)
    # rr is half the difference between k times the k+1'th nearest neighbour distance and
    # the sum of the k nearest neighbour distances
    rr = 0.5 * (k * di[, k+1] - apply(di[, 1:k], MARGIN = 1, FUN = sum))
    # r represents how much further away the (k+1)'th neighbour is than the average distance to the first k neighbours
    r = mean(rr)
    lambda = max(r, 0)
    timer$add('lambda')
    
    #
    # Turns out A and A0 are not used in rest of function so do not execute
    if( FALSE ) {
      # The numerator is the difference between the k-nearest-neighbour differences and the (k+1)'th nearest
      # neighbour difference
      numerator = t(di[, k+1] - t(di))
      # Each row of the denominator is the difference between k times the k+1'th nearest neighbour distance and
      # the sum of the k nearest neighbour distances
      denominator = matrix(rep(2 * rr, k+1), nrow = num)
      # A will be a sparse matrix. The non-zero entries reflect the proportion of the distance to the k nearest
      # neighbours that the i'th nearest neighbour distance is.
      A = array(0, c(num, num))
      # a indexes the rows and columns of A
      a = matrix(rep(1:num, k+1), nrow = num)
      # Index the k+1 nearest neighbours
      id = idx[, 2:(k + 2)]
      # Set the non-zero elements of A
      A[cbind(as.vector(a), as.vector(id))] = as.vector(numerator / denominator)
      # Remove any NaNs
      A[is.nan(A)] = 0
      # A0 is symmetricised version of A
      A0 = (A + t(A)) / 2
    }
    
    #
    # Initialise S, of course this is nothing like as described in the paper
    #
    # Remember distX is the average of the kernel distances
    S = max(distX) - distX
    if( return_intermediaries ) intermediaries$S[1,,] = as.matrix(S)
    # Perform network diffusion
    message("Performing network diffusion.")
    S0 = network.diffusion(S, k)
    # Normalise S0 - this will be used as a starting estimate in the optimisation
    S = dn(S0, 'ave')
    if( return_intermediaries ) intermediaries$Snd[1,,] = as.matrix(S)
    timer$add('diffusion')
    
    #
    # Calculate the Laplacian matrix of the graph represented by the adjacency matrix S
    # https://en.wikipedia.org/wiki/Laplacian_matrix
    D0 = diag(apply(S, MARGIN=2, FUN=sum))
    L0 = D0 - S
    #
    # Spectral decomposition of Laplacian - this relates to the largest eigenvalues of S - I_n talked about in the
    # paper
    eig1_res = eig1(L0, c, 0)
    # F_eig1 is L or at least related to L
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    if( return_intermediaries ) {
      intermediaries$Lvec[1,,] = F_eig1
      intermediaries$Lval[1,] = temp_eig1
    }
    timer$add('update.L')
    
    #
    # Perform the iterative optimisation procedure NITER times
    converge = vector()
    for(iter in 1:NITER) {
      #
      message("Iteration: ", iter)
      #
      # Update S
      #
      # Compute the square of the L2 distances between the columns of the eigenvectors
      # distf is num x num
      distf = as.matrix(dist(F_eig1))^2
      # Construct indexes into the distances
      inda = cbind(rep(1:num, num-1), as.vector(idx[, 2:num]))
      # distX contains the weighted kernel distances
      # Calculate the v_i as in the supplementary material.
      # Note that in the supplementary material, the equation shows that these are subtracted not added
      # The supplementary material is wrong, it should be addition.
      # Also note that the equation uses beta not r in the denominator
      ad = (distX[inda] + lambda * distf[inda]) / 2 / r
      # Convert the vector ad into a matrix
      dim(ad) = c(num, num - 1)
      #
      # call the C function for the optimization
      # This is the piecewise linear and convex optimisation problem mentioned in the supplementary materials.
      # Note that this appears to be parallelisable but is done sequentially in the C code.
      c_input = -t(ad)
      c_output = t(ad)
      ad = t(.Call("projsplx", c_input, c_output))
      #
      # Calculate the adjacency matrix
      A = array(0, c(num, num))
      A[inda] = as.vector(ad)
      # Remove any NaNs
      A[is.nan(A)] = 0
      # Make the adjacency matrix symmetric
      A = (A + t(A)) / 2
      A
      # S is a smoothed version of the old S with the new adjacency
      S = (1 - beta) * S + beta * A
      if( return_intermediaries ) intermediaries$S[iter + 1,,] = as.matrix(S)
      timer$add('update.S')
      # do similarity enhancement by diffusion
      S = network.diffusion(S, k)
      if( return_intermediaries ) intermediaries$Snd[iter + 1,,] = as.matrix(S)
      timer$add('diffusion')
      
      #
      # Update L
      #
      # Note that in the code L is the Laplacian matrix and F_eig1 corresponds to L in the paper
      #
      D = diag(apply(S, MARGIN = 2, FUN = sum))
      L = D - S
      F_old = F_eig1
      eig1_res = eig1(L, c, 0)
      F_eig1 = eig1_res$eigvec
      temp_eig1 = eig1_res$eigval
      ev_eig1 = eig1_res$eigval_full
      if( return_intermediaries ) {
        intermediaries$Lvec[iter + 1,,] = F_eig1
        intermediaries$Lval[iter + 1,] = temp_eig1
      }
      timer$add('update.L')
      
      #
      # Update weights
      #
      DD <- calc.DD(cl, D_Kernels, S)
      alphaK0 = umkl(DD)
      alphaK0 = alphaK0 / sum(alphaK0)
      # Smoothed update of the alphaK parameterised by beta
      alphaK = (1 - beta) * alphaK + beta * alphaK0
      alphaK = alphaK / sum(alphaK)
      if( return_intermediaries ) intermediaries$alphaK[iter + 1,] = alphaK
      timer$add('update.weights')
      
      #
      # Test for convergence
      #
      # This uses the eignengap criterion (between the c'th and (c+1)'th eigenvalues).
      #
      fn1 = sum(ev_eig1[1:c])
      fn2 = sum(ev_eig1[1:(c+1)])
      converge[iter] = fn2 - fn1
      if (iter < 10) {
        # Heuristic to increase lambda and reduce r if there is an eigengap in the first 9 iterations
        if (ev_eig1[length(ev_eig1)] > 0.000001) {
          lambda = 1.5 * lambda
          r = r / 1.01
        }
      } else {
        # If the convergence criterion is getting worse
        if(converge[iter] > converge[iter-1]) {
          # Use the similarity matrix from the last iteration
          S = S_old
          # Advise user if the convergence test warrants it
          if(converge[iter-1] > 0.2) {
            warning('Maybe you should set a larger value of c.')
          }
          # Break out of iteration loop
          break
        }
      }
      # Retain S as S_old in case our next iteration is not good and we want to use it
      S_old = S
      timer$add('convergence')
      
      #
      # Compute the weighted kernel distances
      #
      distX = Reduce("+", lapply(1:length(D_Kernels), function(i) D_Kernels[[i]] * alphaK[i]))
      if( return_intermediaries ) intermediaries$dists[iter + 1,,] = as.matrix(distX)
      #
      # Order the distances
      idx = order.rows(distX, cl = cl)
      timer$add('sort.distances')
    }
    
    #
    # Calculate the Laplacian matrix
    LF = F_eig1
    D = diag(apply(S, MARGIN = 2, FUN = sum))
    L = D - S
    #
    # compute the eigenvalues and eigenvectors of the Laplacian
    eigen_L = eigen(L)
    U = eigen_L$vectors
    D = eigen_L$values
    timer$add('eigen')
    
    #
    # Run t-SNE on the eigenvectors
    do.tsne <- function(k) {
      U_index = seq(ncol(U), (ncol(U) - no.dim + 1))
      tsne(S,k = no.dim,initial_config=U[, U_index])
    }
    if (length(no.dim) == 1) {
      F_last = do.tsne(k = no.dim)
    } else {
      F_last = lapply(no.dim, do.tsne)
    }
    timer$add('t.SNE.eigen')

    
    #
    # Run k-means clustering
    message("Performing Kmeans.")
    y = kmeans(F_last, c, nstart=200)
    timer$add('k.means')
    
    #
    # Run t-SNE on S
    message("Running t-SNE on S.")
    ydata = tsne(S)
    timer$add('t.SNE.S')
    
    #
    # Stop the cluster
    stopCluster(cl)
    
    # create the structure with the results
    results = list()
    results[["y"]] = y
    results[["S"]] = S
    results[["alphaK"]] = alphaK
    results$timings = timer$get_timings()
    if( return_intermediaries ) {
      results$intermediaries = intermediaries
    }
    
    results_all= append(results_all,list(results))
  }
  names(results_all) <- methodname
  matrix = list("Gaussian.Distance" = D_Kernels[[51]], "FZINB.kernel" = FZINB)  
  
  return(list("results" = results_all, "matrix" = matrix, "execution.time" = execution.time))
}
