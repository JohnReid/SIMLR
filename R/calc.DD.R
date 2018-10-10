
calc.DD.parallel <- function(cl, D_Kernels, S)
  parSapply(
    cl = cl,
    X = D_Kernels,
    # This version of FUN is very close numerically to the one below but seems to affect the
    # visualisation plots. That is SIMLR seems numerically unstable when the weights are perturbed
    # slightly
    # FUN = function(DK) sum((.Machine$double.eps + DK) * (S + .Machine$double.eps)) / ncol(DK))
    FUN = function(DK) mean(apply((.Machine$double.eps + DK) * (S + .Machine$double.eps),
                                  MARGIN = 2,
                                  FUN = sum)))

calc.DD.serial <- function(D_Kernels, S)
  sapply(
    D_Kernels,
    function(DK) mean(apply((.Machine$double.eps + DK) * (S + .Machine$double.eps), MARGIN = 2, FUN = sum)))


calc.DD <- function(cl, D_Kernels, S) {
  DD.parallel <- calc.DD.parallel(cl, D_Kernels, S)
  DD.serial <- calc.DD.serial(D_Kernels, S)
  message('DD difference: ', max(abs(DD.parallel - DD.serial)))
  return(DD.parallel)
}
