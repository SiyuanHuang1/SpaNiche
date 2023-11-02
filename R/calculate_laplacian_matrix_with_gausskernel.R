#' calculate_laplacian_matrix_with_gausskernel
#'
#' @param spa_mat spa_mat
#' @param sigma sigma
#'
#' @return
#' @export
#'
#' @examples
calculate_laplacian_matrix_with_gausskernel = function(
    spa_mat,
    sigma = 1
  ){

  dst2=as.matrix(dist(spa_mat)^2)
  # define the adjacency matrix based on a Gaussian kernel
  # sigma is a parameter for the Gaussian kernel you need to specify
  adj_matrix = exp(-1 * dst2 / (2*sigma^2))

  # calculate the degree matrix
  degree_matrix <- diag(rowSums(adj_matrix))

  # finally calculate the Laplacian matrix
  laplacian_matrix <- degree_matrix - adj_matrix

  # return(laplacian_matrix)
  return(list(laplacian = laplacian_matrix,
              degree = degree_matrix,
              adj = adj_matrix))
}
