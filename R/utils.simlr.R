#' @importFrom magrittr %>%
CumulativeTimer <- R6::R6Class("CumulativeTimer",
  public = list(
    last_time = NULL,
    timings = list(),
    initialize = function() {
      self$reset_last_time()
    },
    reset_last_time = function() { self$last_time <- proc.time() },
    add = function(name) {
      # How long since the last call?
      elapsed = proc.time() - self$last_time
      # If we already have an elapsed time for the name, add it
      if( name %in% names(self$timings) ) {
        elapsed = self$timings[[name]] + elapsed
      }
      # Assign to our list
      self$timings[[name]] = elapsed
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
"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {

    # set the needed parameters
    if(is.na(c)) {
        c = dim(A)[1]
    }
    if(c>dim(A)[1]) {
        c = dim(A)[1]
    }
    if(is.na(isMax)) {
        isMax = 1
    }
    if(is.na(isSym)) {
        isSym = 1
    }

    # compute the eigenvalues and eigenvectors of A
    if(isSym==1) {
        eigen_A = eigen(A,symmetric=TRUE)
    }
    else {
        eigen_A = eigen(A)
    }
    v = eigen_A$vectors
    d = eigen_A$values

    # sort the eigenvectors
    if(isMax == 0) {
        eigen_A_sorted = sort(d,index.return=TRUE)
    }
    else {
        eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
    }
    d1 = eigen_A_sorted$x
    idx = eigen_A_sorted$ix
    idx1 = idx[1:c]

    # compute the results
    eigval = d[idx1]
    eigvec = Re(v[,idx1])
    eigval_full = d[idx]

    return(list(eigval=eigval, eigvec=eigvec, eigval_full=eigval_full))

}


#' Compute the L2 distances between the rows of two matrices.
#'
#' Might as well use `as.matrix(dist(a))^2`
#'
#' @param a: First matrix
#' @param b: Second matrix (same shape as a). If NULL will be taken as equivalent to a.
#'
"L2_distance_1" <- function( a, b = NULL ) {
    # Add a row of zeros to a if it only has one row
    if(dim(a)[1] == 1) {
        a = rbind(a, rep(0, dim(a)[2]))
    }
    # Take L2 norm of each column
    aa = apply(a, MARGIN=2, FUN=function(x) sum(x^2))
    # If b is given then use it, otherwise use calculation of a
    if( is.null(b) ) {
        b = a
        bb = aa
    } else {
        # Add a row of zeros to b if it only has one row
        if(dim(b)[1] == 1) {
            b = rbind(b, rep(0, dim(b)[2]))
        }
        # Take L2 norm of each column
        bb = apply(b, MARGIN=2, FUN=function(x) sum(x^2))
    }
    # Calculate distances
    d = Re(outer(aa, bb, FUN = '+') - 2 * crossprod(a, b))
    # Make sure all distances are non-negative
    d[d < 0] = 0
    # Make sure diagonal is exactly 0
    diag(d) = 0
    #
    return(d)
}


#' umkl function
#'
#' TODO: It is unclear what this function does
#'
"umkl" = function( D, beta = NA ) {

    # set some parameters
    if(is.na(beta)) {
        beta = 1 / length(D)
    }
    tol = 1e-4
    u = 20
    logU = log(u)

    # compute Hbeta
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P

    betamin = -Inf
    betamax = Inf
    # evaluate whether the perplexity is within tolerance
    Hdiff = H - logU
    tries = 0
    while (abs(Hdiff) > tol && tries < 30) {
        #if not, increase or decrease precision
        if (Hdiff > 0) {
            betamin = beta
            if(abs(betamax)==Inf) {
                beta = beta * 2
            }
            else {
                beta = (beta + betamax) / 2
            }
        }
        else {
            betamax = beta
            if(abs(betamin)==Inf) {
                beta = beta * 2
            }
            else {
                beta = (beta + betamin) / 2
            }
        }
        # compute the new values
        res_hbeta = Hbeta(D, beta)
        H = res_hbeta$H
        thisP = res_hbeta$P
        Hdiff = H - logU
        tries = tries + 1
    }

    return(thisP)

}


"Hbeta" = function( D, beta ) {

    D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
    P = exp(-D * beta)
    sumP = sum(P)
    H = log(sumP) + beta * sum(D * P) / sumP
    P = P / sumP

    return(list(H=H,P=P))

}


#' Sort each row of X and return the sorted matrix and the ordering vectors
#'
sort.rows = function(X) {
    sorted = array(0, c(nrow(X), ncol(X)))
    idx = array(0, c(nrow(X), ncol(X)))
    for(i in 1:nrow(X)) {
        res = sort(X[i,], index.return = TRUE)
        sorted[i,] = res$x
        idx[i,] = res$ix
    }
    return(list(sorted = sorted, idx = idx))
}
