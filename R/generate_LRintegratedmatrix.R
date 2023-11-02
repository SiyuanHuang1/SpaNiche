#' generate_LRintegratedmatrix
#'
#' @param neighbor_matrix A matrix representing the relationship between neighboring spatial locations.
#' @param neighbor_role A matrix (spot by gene).
#' @param center_role A matrix (spot by gene).
#' @param center_is_ligand A logical variable indicating if the central role represents a ligand. If TRUE, the central role is a ligand; if FALSE, it is a receptor.
#'
#' @return
#' @export
#'
#' @examples
generate_LRintegratedmatrix=function(
    neighbor_matrix,
    neighbor_role,
    center_role,
    center_is_ligand
){
  print("generating L-R integrated matrix:")

  if(is.matrix(neighbor_matrix) &
     is.matrix(neighbor_role) &
     is.matrix(center_role)){
    print("All data is matrix. ok!")
  } else {
    stop("All inputs should be matrices!")
  }

  if (identical(rownames(neighbor_matrix),colnames(neighbor_matrix)) &
      identical(rownames(neighbor_matrix),rownames(neighbor_role)) &
      identical(rownames(neighbor_matrix),rownames(center_role))) {
    print("The order of spatial barcodes is identical. ok!")
  } else {
    stop("The order of spatial barcodes is not identical across matrices!")
  }

  neighbor_of_spot_multiplied_by_neighbor_role_mean=neighbor_matrix %*% neighbor_role
  spot_multiplied_by_center_role=center_role
  region_multiplied_by_LR=neighbor_of_spot_multiplied_by_neighbor_role_mean * spot_multiplied_by_center_role

  if(center_is_ligand){
    LR_name=paste(colnames(center_role),colnames(neighbor_role),sep = "->")
    LR_name=paste0(LR_name," (CisL, NisR)")
  } else {
    LR_name=paste(colnames(center_role),colnames(neighbor_role),sep = "<-")
    LR_name=paste0(LR_name," (CisR, NisL)")
  }

  colnames(region_multiplied_by_LR) = LR_name

  print("Finished.")
  print("------------------------------------")
  return(region_multiplied_by_LR)
}
