#' define_initial_weight
#'
#' @param data.list data.list
#'
#' @return
#' @export
#'
#' @examples
define_initial_weight = function(data.list){

  mean_matrix_norm_v = c()
  element_number_v = c()

  for (i in 1:length(data.list)) {
    one_matrix=data.list[[i]]
    matrix_norm <- norm(one_matrix, "F")

    element_num = dim(one_matrix)[1] * dim(one_matrix)[2]
    element_number_v = append(element_number_v,element_num)

    mean_matrix_norm = matrix_norm / element_num
    mean_matrix_norm_v = append(mean_matrix_norm_v,mean_matrix_norm)
  }

  theta_v = (max(mean_matrix_norm_v) / mean_matrix_norm_v) *
    (max(element_number_v) / element_number_v)
  theta_v
}
