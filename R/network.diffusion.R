#' Perform network diffusion of K steps over the network A
#'
#' This is called similarity enhancement by diffusion in the paper
#'
#' @param A: The distance matrix
#' @param K: The number of nearest neighbours to consider
#' @param scale_by_DD: Scale W by DD (included in standard network.diffusion() but not in network.diffusion.numc()
#'
#' @keywords internal
#'
network.diffusion <- function( A, K, scale_by_DD = TRUE ) {

    # set the values of the diagonal of A to 0
    diag(A) = 0

    # compute the dominate set for A and K
    P = dominate.set(abs(A), min(K, nrow(A) - 1)) * sign(A)

    # sum the absolute value of each row of P
    DD = apply(abs(P), MARGIN = 1, FUN = sum)

    # set the diagonal of P to be DD + 1
    diag(P) = DD + 1

    # compute the transition field of P
    P = transition.fields(P)

    # compute the eigenvalues and eigenvectors of P
    eigen_P = eigen(P)
    U = eigen_P$vectors
    D = eigen_P$values

    # set to d the real part of the diagonal of D
    d = Re(D + .Machine$double.eps)

    # perform the diffusion
    alpha = 0.8
    beta = 2
    d = ((1-alpha)*d)/(1-alpha*d^beta)

    # set D to be a diagonal matrix of the real part of d
    D = diag(Re(d))

    # finally compute W
    W = U %*% D %*% t(U)
    diagonal_matrix = diag(rep(1, nrow(W)))
    W = (W * (1 - diagonal_matrix)) / apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=(1-diag(W))})
    if( scale_by_DD ) {
      # This line is missing in network.diffusion.numc()
      W = diag(DD) %*% W
    }
    # Ensure W symmetric
    W = (W + t(W)) / 2
    # Ensure all W are non-negative
    W[W < 0] = 0

    return(W)
}


#' Compute the dominate set for the matrix aff.matrix and NR.OF.KNN
#'
#' Creates a copy of `aff.matrix` such that every element smaller than the k'th
#' largest in each row is zero'ed. A symmetric version of this matrix is then returned.
#'
#' @keywords internal
#'
dominate.set <- function( aff.matrix, NR.OF.KNN ) {
  PNN.mine <- apply(aff.matrix, MARGIN = 1, FUN = purrr::partial(zero_vec, k = NR.OF.KNN))
  return((t(PNN.mine) + PNN.mine) / 2)
}



#' compute the transition field of the given matrix
#'
#' @keywords internal
#'
transition.fields <- function( W )
{
    # get any index of columns with all 0s
    zero.index = which(apply(W,MARGIN=1,FUN=sum)==0)

    # compute the transition fields
    W = dn(W, 'ave')

    w = sqrt(apply(abs(W),MARGIN=2,FUN=sum)+.Machine$double.eps)
    W = W / t(apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=w}))
    W = W %*% t(W)

    # set to 0 the elements of zero.index
    W[zero.index,] = 0
    W[,zero.index] = 0

    return(W)
}


#' Normalizes a symmetric kernel w
#'
#' @param: w The symmetric kernel
#' @param: type The normalisation type
#'
#' @keywords internal
#'
dn = function( w, type ) {
    #
    # Compute the sums of the columns
    D = apply(w, MARGIN = 2, FUN = sum)
    #
    # type "ave" returns D^-1*W
    if(type=="ave") {
        wn = diag(1 / D) %*% w
    }
    # type "gph" returns D^-1/2*W*D^-1/2
    else if(type=="gph") {
        D_temp = diag(1 / sqrt(D))
        wn = D_temp %*% (w %*% D_temp)
    }
    else {
        stop("Invalid type!")
    }
    return(wn)
}
