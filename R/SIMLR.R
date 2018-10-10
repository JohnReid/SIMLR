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
#' @importFrom parallel parLapply parSapply parRapply parApply
#' @importFrom stats dnorm kmeans pbeta rnorm
#' @importFrom methods is
#' @import Matrix
#' @useDynLib SIMLR projsplx
#'
SIMLR <- function(X,
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
  # Start a cluster for parallel tasks
  cl = start_cluster( cores.ratio )

  #
  # compute the kernel distances
  message("Computing the multiple Kernels.")
  D_Kernels = multiple.kernel(t(X), cl = cl)
  #
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
  # Compute the execution time
  execution.time = proc.time() - ptm

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


#' Apply the SIMLR algorithm to a data set
#'
#' @keywords internal
#'
run_SIMLR <- function(.data)
{
  message("Running SIMLR")
  return(SIMLR(X = .data$in_X, c = .data$n_clust, return_intermediaries = TRUE))
}


#' Summarise the results of the SIMLR algorithm
#
#' @import ggplot2
#'
#' @keywords internal
#'
summarise_SIMLR <- function(
                            res,
                            .data,
                            data_set,
                            output_dir = file.path('output', data_set),
                            max_intermediaries = 6,
                            max_samples = 300)
{
  #
  # Create output directory if needs be
  output_file <- purrr::partial(file.path, output_dir)
  dir.create(output_file('.'), recursive = TRUE, showWarnings = FALSE)

  #
  # Save results to output
  saveRDS(res, output_file('SIMLR-results.rds'))

  #
  # Number of samples
  num <- nrow(res$S)
  sample_idxs <- 1:num
  if( ! is.null(max_samples) && num > max_samples ) {
    # Maintain order as this tends to improve the plots
    sample_idxs <- sort(sample(num, size = max_samples))
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
  # A function to find out how many intermediate iterations we have
  find_last_non_na <- function(x) max(which(! is.na(x)))

  #
  # A function to plot an intermediary similarity/distance matrix
  plot_intermediary <- function(X) {
    last_iter <- find_last_non_na(X[, 1, 1])
    plot_list <- lapply(
                        approx_spaced_integers(1, last_iter, max_intermediaries),
                        function(iter) {
                          similarity.heatmap(X[iter, sample_idxs, sample_idxs],
                                             label = stringr::str_c('label ', .data$true_labs[sample_idxs, 1]),
                                             # cluster = stringr::str_c('cluster ', res$y$cluster[sample_idxs]),
                                             annotation_legend = FALSE,
                                             annotation_names_col = FALSE,
                                             cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             legend = FALSE,
                                             main = iter,
                                             silent = TRUE)[[4]]
                        })
    do.call(gridExtra::grid.arrange, plot_list)
  }

  #
  # Make a grid of the intermediate S
  message('Plotting intermediate S')
  ggsave(output_file("S-intermediaries.pdf"), plot_intermediary(res$intermediaries$S))

  #
  # Make a grid of the intermediate S after network diffusion
  message('Plotting diffused intermediate S')
  ggsave(output_file("S-diffused-intermediaries.pdf"), plot_intermediary(res$intermediaries$Snd))

  #
  # Make a grid of the intermediate distances
  message('Plotting intermediate distances')
  ggsave(output_file("dists-intermediaries.pdf"), plot_intermediary(res$intermediaries$dists))

  #
  # Plot the intermediate weights
  message('Plotting intermediate weights')
  #
  # Get a mapping from kernel indices to parameters
  kernels <- kernel_param_map()
  #
  # Melt the intermediate weights
  iter <- find_last_non_na(res$intermediaries$alphaK[, 1])
  alphaK <-
    reshape2::melt(
                   res$intermediaries$alphaK[1:iter,],
                   varnames=c('iter', 'kernel'),
                   value.name = 'weight') %>%
  dplyr::left_join(kernels)
#
# Make the plot
ggplot(alphaK, aes(x = iter, y = weight, linetype = k, colour = sigma)) + geom_line()
ggsave(output_file("alphaK-intermediaries.pdf"))

#
# Plot the eigenvalues
iter <- find_last_non_na(res$intermediaries$Lval[, 1])
eigvals <-
  reshape2::melt(
                 res$intermediaries$Lval[1:iter,],
                 varnames=c('iter', 'eigenvector'),
                 value.name = 'eigenvalue')
#
# Make the plot
ggplot(eigvals, aes(x = iter, y = eigenvalue, group = eigenvector)) + geom_line()
ggsave(output_file("L-eigenvalues-intermediaries.pdf"))

#
# Show NMI
message('NMI: ', nmi_1)

invisible(res)
}
