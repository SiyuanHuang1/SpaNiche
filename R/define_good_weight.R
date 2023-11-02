#' define_good_weight
#'
#' @param X X
#' @param Wzero Wzero
#' @param Hzero Hzero
#' @param theta_v theta_v
#' @param weight_input weight_input
#'
#' @return
#' @export
#'
#' @examples
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
