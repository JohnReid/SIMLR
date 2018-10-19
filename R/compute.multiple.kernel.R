#' Compute and returns the multiple kernel distances
#'
#' Compute the kernel distances for the hard-wired range sigma (bandwidth scaling) and k (nearest-neighbours)
#'
#' @param x The data (samples x features)
#' @param calc.dists Normalise kernels and convert to kernel distances
#'
#' @keywords internal
#'
multiple.kernel <- function(x, cl, calc.dists = TRUE, offset = 1, dist_power = 2) {
  D_Kernels <- multiple.unnorm.kernels(x, cl = cl, dist_power = dist_power)
  if (calc.dists) {
    D_Kernels <- norm.and.calc.dists(D_Kernels, offset = offset)
  }
  #
  # Return the distances/kernels
  return(D_Kernels)
}


#' Normalise the kernels and convert to kernel distances
#'
#' @param D_Kernels The unnormalised kernel (Gram) matrices
#'
#' @keywords internal
#'
norm.and.calc.dists <- function(D_Kernels, offset = 1) {
  #
  # Normalise the kernels, calculate the kernel distances and convert to sparse matrices
  message("Calculating distances.")
  # for (i in 1:length(D_Kernels)) {
  #     D_Kernels[[i]] = Matrix(kernel.distance.2(kernel.normalise(D_Kernels[[i]])), sparse=TRUE, doDiag=FALSE)
  # }
  D_Kernels <- lapply(
    D_Kernels,
    function(G) kernel.distance.2(kernel.normalise(G, offset = offset))
  )
  return(D_Kernels)
}


#' Compute and returns multiple unnormalised kernels
#'
#' Compute the kernels for the hard-wired range sigma (bandwidth scaling) and k (nearest-neighbours)
#'
#' @param x The data (samples x features)
#' @param cl Cluster to run parallel tasks on
#' @param dist_power The power to raise the distances
#'
#' @keywords internal
#'
multiple.unnorm.kernels <- function(x, cl, dist_power = 2) {
  message("Calculating kernels.")
  #
  # Kernel parameters
  sigma <- default_sigma()
  allk <- default_k()

  #
  # compute and sort Diff
  Diff <- dist2(x)^dist_power # Diff is the square of the squared distance (i.e. power of 4)
  Diff_sort <- t(apply(Diff, MARGIN = 2, FUN = sort)) # Sort the rows to help with kNN later

  # The parallel apply runs over all the k for the kNN
  D_Kernels <- unlist(parLapply(
    cl,
    allk,
    fun = function(k, .Diff_sort = Diff_sort, .Diff = Diff, .sigma = sigma) {
      # Only generate kernels if we have enough data for kNN
      if (k <= ncol(.Diff_sort) - 1) {
        #
        # Calculate the mean of the k-nearest-neighbours,
        # this is mu_i in the paper Eqn. (4)
        TT <- apply(.Diff_sort[, 2:(k + 1)], MARGIN = 1, FUN = mean) + .Machine$double.eps
        #
        # Do an outer average
        # this is (mu_i + mu_j) / 2 in the paper Eqn. (4)
        Sig <- outer(TT, TT, FUN = function(x, y) (x + y) / 2)
        #
        # Ensure every entry is at least machine epsilon
        Sig[Sig < .Machine$double.eps] <- .Machine$double.eps
        #
        # Construct a kernel for each scaling sigma
        sigma_kernels <- lapply(.sigma, FUN = function(sigma) {
          # N.B. .Diff == the squared squared distance (i.e. power of 4)
          return(as(symmpart(dnorm(x = .Diff, mean = 0, sd = sigma * Sig)), 'dspMatrix'))
        })
        return(sigma_kernels)
      }
    }
  ))
  #
  # Return the kernels
  return(D_Kernels)
}


#' Normalise a kernel
#'
#' @keywords internal
#'
kernel.normalise <- function(K, offset = 1) {
  k <- 1 / sqrt(diag(K) + offset)
  return(K * (k %*% t(k)))
}


#' Calculate the kernel distance
#'
#' The kernel distance is the distance in feature space. This function
#' returns half the squared kernel distance for some reason.
#'
#' @param G The kernel (Gram) matrix
#'
#' @keywords internal
#'
kernel.distance.2 <- function(G) {
  # Construct a matrix where each row is the diagonal of the Gram matrix
  G1 <- matrix(rep(diag(G), nrow(G)), nrow = nrow(G))
  # Calculate half the squared distance
  # JR: I'm not sure why the half is there.
  D_Kernels_tmp <- (G1 + t(G1) - 2 * G) / 2
  # Ensure diagonal is 0, it should be anyway
  return(D_Kernels_tmp - diag(diag(D_Kernels_tmp)))
}


#' Compute the squared Euclidean distance
#'
#' @param x1 The data
#' @param x2 The data (assumed to be same as x1 if not given)
#'
#' @keywords internal
#'
dist2 <- function(x1, x2 = NA) {

  # set the parameters for x1
  if (is.na(x2)) {
    x2 <- x1
  }

  # compute the dimension
  n1 <- nrow(x1)
  d1 <- ncol(x1)
  n2 <- nrow(x2)
  d2 <- ncol(x2)
  if (d1 != d2) {
    stop("Data dimension does not match dimension of centres.")
  }

  # compute the distance
  dist <- t(rep(1, n2) %*% t(apply(t(x1^2), MARGIN = 2, FUN = sum))) +
    (rep(1, n1) %*% t(apply(t(x2^2), MARGIN = 2, FUN = sum))) -
    2 * (x1 %*% t(x2))

  return(dist)
}


#' The default parameter set for sigma
#'
#' @keywords internal
#'
default_sigma <- function() seq(2, 1, -0.25)


#' The default parameter set for k
#'
#' @keywords internal
#'
default_k <- function() seq(10, 30, 2)
