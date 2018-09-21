#' Perform network diffusion of K steps over the network A
#'
#' This is called similarity enhancement by diffusion in the paper
#'
#' @param A: The distance matrix
#' @param K: The number of nearest neighbours to consider
#'
"network.diffusion" <- function( A, K ) {

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
    diag(D) = diag(D)[length(diag(D)):1]
    W = diag(DD) %*% W
    # Ensure W wymmetric
    W = (W + t(W)) / 2
    # Ensure all W are non-negative
    W[W < 0] = 0

    return(W)
}


#' Compute the dominate set for the matrix aff.matrix and NR.OF.KNN
"dominate.set" <- function( aff.matrix, NR.OF.KNN ) {

    # create the structure to save the results
    PNN.matrix = array(0,c(nrow(aff.matrix),ncol(aff.matrix)))

    # sort each row of aff.matrix in descending order and saves the sorted 
    # array and a collection of vectors with the original indices
    res.sort = apply(t(aff.matrix),MARGIN=2,FUN=function(x) {return(sort(x, decreasing = TRUE, index.return = TRUE))})
    sorted.aff.matrix = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$x) }))
    sorted.indices = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$ix) }))

    # get the first NR.OF.KNN columns of the sorted array
    res = sorted.aff.matrix[,1:NR.OF.KNN]

    # create a matrix of NR.OF.KNN columns by binding vectors of 
    # integers from 1 to the number of rows/columns of aff.matrix
    inds = array(0,c(nrow(aff.matrix),NR.OF.KNN))
    inds = apply(inds,MARGIN=2,FUN=function(x) {x=1:nrow(aff.matrix)})

    # get the first NR.OF.KNN columns of the indices of aff.matrix
    loc = sorted.indices[,1:NR.OF.KNN]

    # assign to PNN.matrix the sorted indices
    PNN.matrix[(as.vector(loc)-1)*nrow(aff.matrix)+as.vector(inds)] = as.vector(res)

    # compute the final results and return them
    PNN.matrix = (PNN.matrix + t(PNN.matrix))/2

    return(PNN.matrix)

}

# compute the transition field of the given matrix
"transition.fields" <- function( W ) {

    # get any index of columns with all 0s
    zero.index = which(apply(W,MARGIN=1,FUN=sum)==0)

    # compute the transition fields
    W = dn(W,'ave')

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
"dn" = function( w, type ) {
    #
    # Compute the sums of the columns
    D = apply(w, MARGIN=2, FUN=sum)
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
