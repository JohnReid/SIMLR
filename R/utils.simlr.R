#' A class to measure cumulative time across sub-tasks.
#'
#' @importFrom magrittr %>%
#'
#' @keywords internal
#'
CumulativeTimer <-
  R6::R6Class(
    "CumulativeTimer",
    public = list(
      last_time = NULL,
      timings = list(),
      initialize = function() {
        self$reset_last_time()
      },
      reset_last_time = function() {
        self$last_time <- proc.time()
      },
      add = function(name) {
        # How long since the last call?
        elapsed <- proc.time() - self$last_time
        # If we already have an elapsed time for the name, add it
        if (name %in% names(self$timings)) {
          elapsed <- self$timings[[name]] + elapsed
        }
        # Assign to our list
        self$timings[[name]] <- elapsed
        # Update the last time
        self$reset_last_time()
      },
      get_timings = function() {
        as.data.frame(t(as.data.frame(lapply(self$timings, data.matrix)))) %>%
          dplyr::mutate(task = names(self$timings)) %>%
          dplyr::arrange(-elapsed)
      }
    )
  )


#' Compute the eigenvalues and eigenvectors of A
#'
#' @param A The matrix to compute the eigenvalues and vectors of
#' @param c The number of eigenvectors to compute
#' @param isMax Either return the largest or the smallest eigenvalues (and associated vectors)
#' @param isSym Is A symmetric?
#'
#' @keywords internal
#'
eig1 <- function(A, c = NA, isMax = NA, isSym = NA) {

  # set the needed parameters
  if (is.na(c)) {
    c <- dim(A)[1]
  }
  if (c > dim(A)[1]) {
    c <- dim(A)[1]
  }
  if (is.na(isMax)) {
    isMax <- 1
  }
  if (is.na(isSym)) {
    isSym <- 1
  }

  # compute the eigenvalues and eigenvectors of A
  if (isSym == 1) {
    eigen_A <- eigen(A, symmetric = TRUE)
  }
  else {
    eigen_A <- eigen(A)
  }
  v <- eigen_A$vectors
  d <- eigen_A$values

  # sort the eigenvectors
  if (isMax == 0) {
    eigen_A_sorted <- sort(d, index.return = TRUE)
  }
  else {
    eigen_A_sorted <- sort(d, decreasing = TRUE, index.return = TRUE)
  }
  d1 <- eigen_A_sorted$x
  idx <- eigen_A_sorted$ix
  idx1 <- idx[1:c]

  # compute the results
  eigval <- d[idx1]
  eigvec <- Re(v[, idx1])
  eigval_full <- d[idx]

  return(list(eigval = eigval, eigvec = eigvec, eigval_full = eigval_full))
}


#' Compute the L2 distances between the rows of two matrices.
#'
#' Might as well use `as.matrix(dist(a))^2`
#'
#' @param a: First matrix
#' @param b: Second matrix (same shape as a). If NULL will be taken as equivalent to a.
#'
#' @keywords internal
#'
L2_distance_1 <- function(a, b = NULL) {
  # Add a row of zeros to a if it only has one row
  if (dim(a)[1] == 1) {
    a <- rbind(a, rep(0, dim(a)[2]))
  }
  # Take L2 norm of each column
  aa <- apply(a, MARGIN = 2, FUN = function(x) sum(x^2))
  # If b is given then use it, otherwise use calculation of a
  if (is.null(b)) {
    b <- a
    bb <- aa
  } else {
    # Add a row of zeros to b if it only has one row
    if (dim(b)[1] == 1) {
      b <- rbind(b, rep(0, dim(b)[2]))
    }
    # Take L2 norm of each column
    bb <- apply(b, MARGIN = 2, FUN = function(x) sum(x^2))
  }
  # Calculate distances
  d <- Re(outer(aa, bb, FUN = "+") - 2 * crossprod(a, b))
  # Make sure all distances are non-negative
  d[d < 0] <- 0
  # Make sure diagonal is exactly 0
  diag(d) <- 0
  #
  return(d)
}


#' umkl function
#'
#' TODO: It is unclear what this function does
#'
#' @keywords internal
#'
umkl <- function(D, beta = NA) {

  # set some parameters
  if (is.na(beta)) {
    beta <- 1 / length(D)
  }
  tol <- 1e-4
  u <- 20
  logU <- log(u)

  # compute Hbeta
  res_hbeta <- Hbeta(D, beta)
  H <- res_hbeta$H
  thisP <- res_hbeta$P

  betamin <- -Inf
  betamax <- Inf
  # evaluate whether the perplexity is within tolerance
  Hdiff <- H - logU
  tries <- 0
  while (abs(Hdiff) > tol && tries < 30) {
    # if not, increase or decrease precision
    if (Hdiff > 0) {
      betamin <- beta
      if (abs(betamax) == Inf) {
        beta <- beta * 2
      }
      else {
        beta <- (beta + betamax) / 2
      }
    }
    else {
      betamax <- beta
      if (abs(betamin) == Inf) {
        beta <- beta * 2
      }
      else {
        beta <- (beta + betamin) / 2
      }
    }
    # compute the new values
    res_hbeta <- Hbeta(D, beta)
    H <- res_hbeta$H
    thisP <- res_hbeta$P
    Hdiff <- H - logU
    tries <- tries + 1
  }

  return(thisP)
}


Hbeta <- function(D, beta) {
  D <- (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
  P <- exp(-D * beta)
  sumP <- sum(P)
  H <- log(sumP) + beta * sum(D * P) / sumP
  P <- P / sumP

  return(list(H = H, P = P))
}


#' Sort each row of X and return the sorted matrix and the ordering vectors
#'
#' @keywords internal
#'
sort.rows <- function(X, cl, decreasing = FALSE) {
  sorted <- array(0, c(nrow(X), ncol(X)))
  idx <- array(0, c(nrow(X), ncol(X)))
  for (i in 1:nrow(X)) {
    res <- sort(X[i, ], index.return = TRUE, decreasing = decreasing)
    sorted[i, ] <- res$x
    idx[i, ] <- res$ix
  }
  return(list(sorted = sorted, idx = idx))
}

#' Return a vector of the same length as x, where all elements that are less
#' than the k'th largest of x are 0. Other elements are untouched.
#'
#' @keywords internal
#'
zero_vec <- function(x, k) {
  res <- rep(0, length(x))
  partial <- length(x) - k + 1
  idxs <- x >= sort(x, partial = partial)[partial]
  res[idxs] <- x[idxs]
  return(res)
}


#' Order each row of X and return the ordering vectors
#'
#' @keywords internal
#'
order.rows <- function(X, cl) t(parApply(cl, X, MARGIN = 1, FUN = order))


#' Return a vector of integers from from to to that are approximately evenly spaced
#'
#' @keywords internal
#'
approx_spaced_integers <- function(from, to, max_length) {
  result <- from:to
  if (!is.null(max_length) && length(result) > max_length) {
    result <- round(seq(from, to, length.out = max_length))
  }
  return(result)
}


#' Return a data frame containing the kernel parameters.
#'
#' @keywords internal
#'
kernel_param_map <- function() {
  #
  # Hard-code these as they are hard-coded elsewhere
  allk <- seq(10, 30, 2)
  sigma <- seq(2, 1, -0.25)
  return(data.frame(
    kernel = 1:55,
    k = factor(rep(allk, each = length(sigma))),
    sigma = factor(rep(sigma, length(allk)))
  ))
}


#' Calculate how many cores from the ratio
#'
#' @keywords internal
#'
cores_from_ratio <- function(cores.ratio) {
  cores <- as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores <- 1
  }
  return(cores)
}


#' Set up the cluster to parallelise the kernel calculations
#'
#' @param cores.ratio Proportional of all possible cores - 1 to use.
#'
#' @keywords internal
#'
start_cluster <- function(num.cores = 0, cores.ratio = 1) {
  if (!num.cores) {
    # Choose how many cores we will use for parallelisation
    num.cores <- cores_from_ratio(cores.ratio)
  }
  cl <- makeCluster(num.cores)
  clusterEvalQ(cl, {
    library(Matrix)
  })
  return(cl)
}


#' Sum the columns
#'
#' @keywords internal
#'
col_sums <- function(W, method = "apply")
  switch(method,
    "apply" = apply(W, MARGIN = 2, FUN = sum),
    "colSums" = colSums(W),
    stop("Unknown method")
  )


#' Sum the rows
#'
#' @keywords internal
#'
row_sums <- function(W, method = "apply")
  switch(method,
    "apply" = apply(W, MARGIN = 1, FUN = sum),
    "rowSums" = rowSums(W),
    stop("Unknown method")
  )


#' Scale the columns by w
#'
#' @keywords internal
#'
scale_cols <- function(W, w, method = "denom")
  switch(method,
    orig = W / t(apply(array(0, c(nrow(W), ncol(W))), MARGIN = 2, FUN = function(x) {
      x <- w
    })),
    denom = {
      denom <- t(matrix(rep(w, ncol(W)), ncol = nrow(W)))
      W / denom
    },
    dense = as.matrix(W) %*% Diagonal(x = 1 / w),
    sparse = W %*% Diagonal(x = 1 / w),
    stop("Unknown method")
  )


#' Scale the rows by w.
#'
#' @keywords internal
#'
scale_rows <- function(W, w, method = "sparse")
  switch(method,
    dense = diag(1 / w) %*% W,
    sparse = Diagonal(x = 1 / w) %*% W,
    stop("Unknown method")
  )


#' Calculate eigenvalues of a symmetric matrix P.
#'
#' @importFrom Rspectra eigs
#'
#' @keywords internal
#'
calc_eigs <- function(P, method = "eigen")
  switch(method,
    eigen = eigen(P, symmetric = TRUE),
    Rspectra = eigs_sym(as(P, "dgCMatrix"), ncol(P) - 1), # Will use eigen if asked for all eigendimensions
    stop("Unknown method")
  )


#' Divide rows by diagonal
#'
#' @keywords internal
#'
divide_rows_by_diag <- function(W, method = "sweep") {
  diagW <- diag(W)
  diag(W) <- 0
  switch(method,
    sweep = sweep(W, MARGIN = 1, STATS = 1 - diagW, FUN = '/'),
    apply = W / apply(array(0, c(nrow(W), ncol(W))), MARGIN = 2, FUN = function(x) x <- (1 - diagW)),
    stop("Unknown method")
  )
}


#' Multiply rows, i.e. do \code{diag(v) %*% X}
#'
#' @keywords internal
#'
multiply_rows <- function(v, X, method = "sweep")
  switch(method,
    sweep = sweep(X, MARGIN = 1, STATS = v, FUN = '*'),
    diag = diag(v) %*% X,
    stop("Unknown method")
  )


#' Multiply cols, i.e. do \code{X %*% diag(v)}
#'
#' @keywords internal
#'
multiply_cols <- function(v, X, method = "sweep")
  switch(method,
    sweep = sweep(X, MARGIN = 2, STATS = v, FUN = '*'),
    diag = X %*% diag(v),
    stop("Unknown method")
  )
