#' perform the SIMLR clustering algorithm
#'
#' @title SIMLR
#'
#' @examples
#' SIMLR(X = BuettnerFlorian$in_X, c = BuettnerFlorian$n_clust, cores.ratio = 0)
#'
#' @param X an (m x n) data matrix of gene expression measurements of individual cells or
#'   an object of class SCESet
#' @param c number of clusters to be estimated over X
#' @param no.dim number of dimensions
#' @param k tuning parameter
#' @param if.impute should I traspose the input data?
#' @param normalize should I normalize the input data?
#' @param cores.ratio ratio of the number of cores to be used when computing the multi-kernel
#' @param return_intermediaries Return intermediate values of S
#'
#' @return clusters the cells based on SIMLR and their similarities
#'
#' @return list of 8 elements describing the clusters obtained by SIMLR, of which y are the resulting clusters:
#'    y = results of k-means clusterings,
#'    S = similarities computed by SIMLR,
#'    F = results from network diffiusion,
#'    ydata = data referring the the results by k-means,
#'    alphaK = clustering coefficients,
#'    execution.time = execution time of the present run,
#'    timings = execution times of sub-tasks,
#'    converge = iterative convergence values by T-SNE,
#'    LF = parameters of the clustering
#'
#' @export SIMLR
#' @importFrom parallel stopCluster makeCluster detectCores clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom stats dnorm kmeans pbeta rnorm
#' @importFrom methods is
#' @import Matrix
#' @useDynLib SIMLR projsplx
#'
"SIMLR" <- function(X,
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
        cat("X is an SCESet, converting to input matrix.\n")
        X = X@assayData$exprs
    }

    # set some parameters
    NITER = 30
    num = ncol(X)
    beta = 0.8

    # Measure execution times of sub-tasks
    timer = CumulativeTimer$new()

    # set any required parameter to the defaults
    if(is.na(no.dim)) {
        no.dim = c
    }

    # Check the return intermediaries parameter
    if( return_intermediaries ) {
      intermediaries = list(
        S = array(NA, dim = c(NITER+1, num, num)),
        alphaK = array(NA, dim = c(NITER, 11 * 5)),
        dists = array(NA, dim = c(NITER+1, num, num)))
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

    # Remember the time to calculate execution time later
    ptm = proc.time()

    #
    # compute the kernel distances
    cat("Computing the multiple Kernels.\n")
    D_Kernels = multiple.kernel(t(X), cores.ratio)
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
    temp = sort.rows(distX)
    distX1 = temp$sorted
    idx = temp$idx
    timer$add('sort.distances')

    #
    # Use data to determine lambda
    #
    # di contains the distances to the (k+1) nearest neighbours
    di = distX1[,2:(k+2)]
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
    # Perform network diffusion
    cat("Performing network diffusion.\n")
    # Remember distX is the average of the kernel distances
    S0 = network.diffusion(max(distX) - distX, k)
    # Normalise S0 - this will be used as a starting estimate in the optimisation
    S = dn(S0, 'ave')
    if( return_intermediaries ) intermediaries$S[1,,] = S
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
    evs_eig1 = eig1_res$eigval_full
    timer$add('update.L')

    #
    # Perform the iterative optimisation procedure NITER times
    converge = vector()
    for(iter in 1:NITER) {

        cat("Iteration: ", iter, "\n")

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
        timer$add('update.S')
        # do similarity enhancement by diffusion
        S = network.diffusion(S, k)
        if( return_intermediaries ) intermediaries$S[iter + 1,,] = S
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
        evs_eig1 = cbind(evs_eig1, ev_eig1)
        timer$add('update.L')

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
        # We need the sorted distances updated according to the new weights
        # sort each row of distX into distX1 and retain the ordering vectors in idx
        temp = sort.rows(distX)
        distX1 = temp$sorted
        idx = temp$idx
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
    # Compute the execution time
    execution.time = proc.time() - ptm

    #
    # Run k-means clustering
    cat("Performing Kmeans.\n")
    y = kmeans(F_last, c, nstart=200)
    timer$add('k.means')

    #
    # Run t-SNE on S
    cat("Running t-SNE on S.\n")
    ydata = tsne(S)
    timer$add('t.SNE.S')

    # create the structure with the results
    result = list(
        y = y,
        S = S,
        F = F_last,
        ydata = ydata,
        alphaK = alphaK,
        execution.time = execution.time,
        timings = timer$get_timings(),
        iter = iter,
        converge = converge,
        LF = LF)
    if( return_intermediaries ) {
      result$intermediaries = intermediaries
    }
    return(result)
}


#' Apply the SIMLR algorithm to a data set and summarise the results
#'
apply_SIMLR <- function(
  .data,
  data_set,
  res = NULL,  # Pre-generated results
  output_dir = file.path('output', data_set),
  max_intermediaries = 6)
{
  #
  # run SIMLR
  if( is.null(res) ) {
    message("Running SIMLR")
    res = SIMLR(X = .data$in_X, c = .data$n_clust, return_intermediaries = TRUE)
  }

  #
  # Calculate NMI
  nmi_1 = igraph::compare(.data$true_labs[,1], res$y$cluster, method = "nmi")

  #
  # Report results
  print('Iterations:')
  print(res$iter)
  print(res$execution.time)
  print('Convergence:')
  print(res$converge)
  print('Weights:')
  print(res$alphaK)

  #
  # Scatter plot of dimensionality reduction
  output_file <- purrr::partial(file.path, output_dir)
  dir.create(output_file('.'), recursive = TRUE, showWarnings = FALSE)
  pdf(output_file('scatter.pdf'), width=9, height=6, paper='special')
  plot(res$ydata,
      col = get_palette(3)[.data$true_labs[,1]],
      xlab = "SIMLR component 1",
      ylab = "SIMLR component 2",
      pch = 20,
      main = "SIMILR 2D visualization")
  dev.off()


  #
  # Make a heatmap of S
  # devtools::load_all('../..')
  pdf(output_file('heatmap.pdf'), width=8, height=8, paper='special')
  similarity.heatmap(res$S,
                    label = stringr::str_c('label ', .data$true_labs[,1]),
                    cluster = stringr::str_c('cluster ', res$y$cluster))
  dev.off()

  #
  # Show and save timings
  print(res$timings)
  readr::write_csv(res$timings %>% dplyr::mutate(data.set = data_set, niter = res$iter), output_file('timings.csv'))

  #
  # Make a grid of the intermediate S
  plot_list = list()
  for (iter in approx_spaced_integers(1, res$iter + 1, max_intermediaries)) {
    ph <- similarity.heatmap(res$intermediaries$S[iter,,],
                             label = stringr::str_c('label ', .data$true_labs[,1]),
                             # cluster = stringr::str_c('cluster ', res$y$cluster),
                             annotation_legend = FALSE,
                             annotation_names_col = FALSE,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             legend = FALSE,
                             main = iter)
    plot_list[[length(plot_list) + 1]] <- ph[[4]]  # to save each plot into a list. note the [[4]]
  }
  ggplot2::ggsave(output_file("S-intermediaries.pdf"), do.call(gridExtra::grid.arrange, plot_list))

  #
  # Make a grid of the intermediate distances
  plot_list = list()
  for (iter in approx_spaced_integers(1, res$iter, max_intermediaries)) {
    print(iter)
    ph <- similarity.heatmap(res$intermediaries$dists[iter,,],
                             label = stringr::str_c('label ', .data$true_labs[,1]),
                             # cluster = stringr::str_c('cluster ', res$y$cluster),
                             annotation_legend = FALSE,
                             annotation_names_col = FALSE,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             legend = FALSE,
                             main = iter)
    plot_list[[length(plot_list) + 1]] <- ph[[4]]  # to save each plot into a list. note the [[4]]
  }
  ggplot2::ggsave(output_file("dists-intermediaries.pdf"), do.call(gridExtra::grid.arrange, plot_list))

  #
  # Plot the intermediate weights
  #
  # Hard-code these as they are hard-coded elsewhere
  allk <- seq(10, 30, 2)
  sigma <- seq(2, 1, -0.25)
  kernels <- data.frame(
    kernel = 1:55,
    k = factor(rep(allk, each = length(sigma))),
    sigma = factor(rep(sigma, length(allk))))
  #
  # Melt the intermediate weights
  alphaK <-
    reshape2::melt(
      res$intermediaries$alphaK[1:(res$iter+1),],
      varnames=c('iter', 'kernel'),
      value.name = 'weight') %>%
    left_join(kernels)
  alphaK %>% sample_n(15)
  #
  # Make the plot
  ggplot(alphaK, aes(x = iter, y = weight, linetype = k, colour = sigma)) + geom_line()
  ggplot2::ggsave(output_file("alphaK-intermediaries.pdf"))

  #
  # Show NMI
  message('NMI: ', nmi_1)

  return(res)
}
