#' Perform network diffusion using K neareast neighbours over the network A
#'
#' This is called similarity enhancement by diffusion in the paper
#'
#' @param A: The distance matrix
#' @param K: The number of nearest neighbours to consider
#' @param scale_by_DD: Scale W by DD (included in standard network.diffusion() but not in network.diffusion.numc()
#'
#' @keywords internal
#'
network.diffusion <- function(A, K, scale_by_DD = TRUE) {

  # Check our assumptions about arguments are correct
  stopifnot(all(A >= 0))

  # set the values of the diagonal of A to 0
  diag(A) <- 0

  # compute the dominate set for A and K
  P <- dominate.set(abs(A), min(K, nrow(A) - 1)) * sign(A)

  # sum the absolute value of each row of P
  DD <- row_sums(abs(P), method = "apply")

  # set the diagonal of P to be DD + 1
  diag(P) <- DD + 1

  # compute the transition field of P
  P <- transition.fields(P)

  # compute the eigenvalues and eigenvectors of P
  eigen_P <- calc_eigs(P, method = "eigen")
  U <- eigen_P$vectors
  D <- eigen_P$values

  # set to d the real part of the diagonal of D + some eps
  d <- Re(D + .Machine$double.eps)

  # Adjustment: down-weight the smaller eigenvalues. Is this the diffusion?
  alpha <- 0.8
  beta <- 2
  d <- ((1 - alpha) * d) / (1 - alpha * d^beta)

  # Reconstruct W from the (adjusted) eigendecomposition
  W <- U %*% diag(d) %*% t(U)

  # Do some (weird?) calculation that zeros the diagonal and divides
  # every other element by the diagonal element on its row
  W <- divide_rows_by_diag(W)

  # Scale rows by DD if requested
  # This is missing in network.diffusion.numc()
  if (scale_by_DD) {
    W <- multiply_rows(DD, W)
  }

  # Ensure W symmetric
  W <- (W + t(W)) / 2

  # Ensure all W are non-negative
  W[W < 0] <- 0

  return(W)
}


#' Calculate the nearest neighbours of each sample.
#'
#' Creates a copy of `aff.matrix` such that every element smaller than the k'th
#' largest in each row is zero'ed. A symmetric version of this matrix is then returned.
#'
#' @keywords internal
#'
dominate.set <- function(aff.matrix, NR.OF.KNN) {
  PNN <- apply(aff.matrix, MARGIN = 1, FUN = purrr::partial(zero_vec, k = NR.OF.KNN))
  return(as((t(PNN) + PNN) / 2, "dsCMatrix"))
}


#' Compute the transition field of the given matrix
#'
#' @importFrom Matrix t tcrossprod Diagonal colSums
#'
#' @keywords internal
#'
transition.fields <- function(W) {
  # the indexes of rows that sum to 0
  # I don't believe any rows can sum to 1 as each row contains
  # all elements of S that are not smaller than the row's k'th
  # element
  zero.index <- which(apply(W, MARGIN = 1, FUN = sum) == 0)
  #
  # Double check our assumption is correct....
  stopifnot(length(zero.index) == 0)
  #
  # Normalise W
  W <- dn(W, "ave")
  #
  # Divide each element by the square root of the sum of the absolute value of its column
  w <- sqrt(col_sums(abs(W), "apply") + .Machine$double.eps)
  W <- scale_cols(W, w)
  #
  # Cross product W = W %*% t(W)
  W <- tcrossprod(W)
  #
  # set the elements of zero.index to 0
  if (length(zero.index)) {
    W[zero.index, ] <- 0
    W[, zero.index] <- 0
  }
  #
  return(W)
}


#' Normalizes a symmetric kernel w.
#'
#' @param: w The symmetric kernel
#' @param: type The normalisation type
#'    \enumerate{
#'      \item 'ave' Scales the rows by the column sums
#'      \item 'gph' Scales each element by the square root of its column and row sums
#'    }
#'
#' @keywords internal
#'
dn <- function(W, type) {
  #
  # Scaling factors
  w <- col_sums(W, "apply")
  #
  # Normalisation depending on type.
  switch(type,
    ave = scale_rows(W, w),
    gph = {
      w_sqrt <- Diagonal(x = 1 / sqrt(w))
      w_sqrt %*% (W %*% w_sqrt)
    },
    stop("Invalid normalisation type!")
  )
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
