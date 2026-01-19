#' define_good_weight: define adaptive weights based on reconstruction error
#'
#' @description
#' This function computes an adaptive weight for each matrix based on
#' its reconstruction error under a given matrix factorization.
#' @param X A list of numeric matrices. Each matrix represents an observed data matrix to be factorized.
#' @param Wzero Wzero
#' @param Hzero Hzero
#' @param theta_v A numeric vector of initial weights for each matrix.
#' @param weight_input A numeric vector specifying the target or prior weights for each matrix.
#'
#' @return A numeric vector of adaptive weights, with length equal to length(X).
#' @keywords internal
define_good_weight = function(X,Wzero,Hzero,theta_v,weight_input){

  rho_v = c()

  for (i in 1:length(X)) {
    matrix_norm <- norm(X[[i]] - Wzero %*% Hzero[[i]], "F")
    matrix_norm_square = matrix_norm^2
    one_delta = weight_input[i] / (matrix_norm_square*theta_v[i])
    one_rho = one_delta * theta_v[i]
    rho_v = append(rho_v,one_rho)
  }

  rho_v
}
