#' The SIMLR clustering algorithm
#'
#' @param X genes x samples expression matrix
#' @param c Number of clusters (parameter name clashes with c()!)
#' @param no.dim Dimensions to reduce to (using t-SNE)
#' @param k Used as a k-nearest-neighbour parameter
#' @return A list:
#'   \describe{
#'     \item{execution.time}{Time taken for core of algorithm}
#'     \item{convergence}{Objective function values used to evaluate convergence}
#'     \item{alphaK}{Weights for weighted sum of kernels}
#'     \item{y}{Result of running k-means}
#'     \item{S}{Learned similarity matrix}
#'     \item{F}{Output from t-SNE application to S with different numbers of dimensions}
#'     \item{LF}{The latent matrix L}
#'     \item{ydata}{Output from t-SNE application to S}
#'   }
#'
"SIMLR" <- function( X, c, no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {

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

    # check the normalize parameter
    if(normalize == TRUE) {
        X = t(X)
        X = X - min(as.vector(X))
        X = X / max(as.vector(X))
        C_mean = as.vector(colMeans(X))
        X = apply(X,MARGIN=1,FUN=function(x) return(x-C_mean))
    }

    # Measure execution times of sub-tasks
    timings = list()
    add.timings = function(name) {
      now = proc.time()
      elapsed = now - last.time
      if( name %in% names(timings) ) {
        elapsed = timings[[name]] + elapsed
      }
      timings[[name]] = elapsed
      assign("timings", value = timings, envir = parent.frame())
      return(now)
    }
    last.time = ptm = proc.time()


    # set some parameters
    NITER = 30
    num = ncol(X)
    r = -1
    beta = 0.8

    # compute the kernel distances
    cat("Computing the multiple Kernels.\n")
    D_Kernels = multiple.kernel(t(X), cores.ratio)
    last.time = add.timings('distances')

    #
    # set up some parameters
    #
    # alphaK looks like a Dirichlet prior: 1 / # categories
    alphaK = 1 / rep(length(D_Kernels), length(D_Kernels))
    #
    # distX is the average of the distances
    distX = Reduce("+", D_Kernels) / length(D_Kernels)
    #
    # sort each row of distX into distX1 and retain the ordering vectors in idx
    temp = sort.rows(distX)
    distX1 = temp$sorted
    idx = temp$idx

    #
    # Use data to determine lambda
    #
    # di contains the distances to the (k+1) nearest neighbours
    di = distX1[,2:(k+2)]
    # rr is half the difference between k times the k+1'th nearest neighbour distance and
    # the sum of the k nearest neighbour distances
    rr = 0.5 * (k * di[, k+1] - apply(di[, 1:k], MARGIN = 1, FUN = sum))
    # r is hard-coded to -1 so is always less than 0
    if(r <= 0) {
        r = mean(rr)
    }
    # r represents how much further away the (k+1)'th neighbour is than the average distance to the first k neighbours
    lambda = max(mean(rr), 0)
    last.time = add.timings('lambda')

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
    # Perform network diffusion
    cat("Performing network diffusion.\n")
    # Remember distX is the average of the kernel distances
    S0 = network.diffusion(max(distX) - distX, k)
    # Normalise S0 - this will be used as a starting estimate in the optimisation
    S = dn(S0, 'ave')
    last.time = add.timings('diffusion')

    #
    # Calculate the Laplacian matrix of the graph represented by the adjacency matrix S
    # https://en.wikipedia.org/wiki/Laplacian_matrix
    D0 = diag(apply(S, MARGIN=2, FUN=sum))
    L0 = D0 - S
    #
    # Spectral decomposition of Laplacian
    eig1_res = eig1(L0, c, 0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    evs_eig1 = eig1_res$eigval_full
    last.time = add.timings('spectral')

    #
    # Perform the iterative optimisation procedure NITER times
    converge = vector()
    for(iter in 1:NITER) {

        cat("Iteration: ", iter, "\n")

        #
        # Update S
        #
        # Compute the L2 distances between the transpose of the eigenvectors
        distf = L2_distance_1(t(F_eig1), t(F_eig1))
        # b contains the indexes (sorting order) vectors for all the same cell
        b = idx[, 2:num]
        inda = cbind(rep(1:num, num-1), as.vector(b))
        ad = (distX[inda] + lambda * distf[inda]) / 2 / r
        dim(ad) = c(num, ncol(b))
        #
        # call the C function for the optimization
        c_input = -t(ad)
        c_output = t(ad)
        ad = t(.Call("projsplx_R", c_input, c_output))
        last.time = add.timings('update.S')
        #
        # calculate the adjacency matrix
        A = array(0, c(num, num))
        A[inda] = as.vector(ad)
        # Remove any NaNs
        A[is.nan(A)] = 0
        # Make the adjacency matrix symmetric
        A = (A + t(A)) / 2
        # S is a smoothed version of the old S with the new adjacency
        S = (1 - beta) * S + beta * A
        # do network diffusion again (this is not mentioned in the paper or supplementary materials
        S = network.diffusion(S, k)
        last.time = add.timings('diffusion')

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
        evs_eig1 = cbind(evs_eig1, ev_eig1)
        last.time = add.timings('spectral')

        #
        # Update weights
        #
        DD = sapply(
            D_Kernels,
            function(DK) mean(apply((.Machine$double.eps + DK) * (S + .Machine$double.eps), MARGIN = 2, FUN = sum)))
        alphaK0 = umkl(DD)
        alphaK0 = alphaK0 / sum(alphaK0)
        # Smoothed update of the alphaK parameterised by beta
        alphaK = (1 - beta) * alphaK + beta * alphaK0
        alphaK = alphaK / sum(alphaK)
        last.time = add.timings('update.weights')

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
        last.time = add.timings('convergence')

        #
        # Compute Kbeta, the weighted kernel distances
        #
        distX = Reduce("+", lapply(1:length(D_Kernels), function(i) D_Kernels[[i]] * alphaK[i]))

        #
        # We need the sorted distances updated according to the new weights
        # sort each row of distX into distX1 and retain the ordering vectors in idx
        temp = sort.rows(distX)
        distX1 = temp$sorted
        idx = temp$idx
        last.time = add.timings('sort.distances')
    }

    #
    # Calculate the Laplacian matrix
    LF = F_eig1
    D = diag(apply(S, MARGIN = 2, FUN = sum))
    L = D - S

    #
    # compute the eigenvalues and eigenvectors of P
    eigen_L = eigen(L)
    U = eigen_L$vectors
    D = eigen_L$values
    last.time = add.timings('eigen')

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
    last.time = add.timings('t.SNE.eigen')

    #
    # Compute the execution time
    execution.time = proc.time() - ptm

    #
    # Run k-means clustering
    cat("Performing Kmeans.\n")
    y = kmeans(F_last, c, nstart=200)
    last.time = add.timings('k.means')

    #
    # Run t-SNE on S
    cat("Running t-SNE on S.\n")
    ydata = tsne(S)
    last.time = add.timings('t.SNE.S')

    # create the structure with the results
    return(list(
        y = y,
        S = S,
        F = F_last,
        ydata = ydata,
        alphaK = alphaK,
        execution.time = execution.time,
        timings = timings,
        converge = converge,
        LF = LF))
}
