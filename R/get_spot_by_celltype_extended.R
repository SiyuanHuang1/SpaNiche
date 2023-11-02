#' get_spot_by_celltype_extended
#'
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". The row names represent the barcode of the spot.
#' @param spot_by_celltype a data frame or matrix representing the abundance of cell types, with spots as rows and cell types as columns. Its row names should be consistent with the row names of 'spatialdf'.
#' @param smoothing_type This parameter has three options: c('loop0'), c('loop0', 'loop1'), and c('loop0', 'loop1', 'loop1+loop2'). 'loop0', 'loop1', and 'loop1+loop2' indicate the range within which data smoothing is performed, denoting: no smoothing, smoothing using one adjacent spot, and smoothing using two surrounding spots, respectively. As it involves the merging of matrices, the returned matrix's dimension varies depending on the chosen option. The dimensions corresponding to the three options are 1x, 2x, and 3x of the original dimensions, respectively.
#' @param distance_thre This parameter represents distance, defining the central distance from loop0 to loop1, and from loop0 to loop2. The length of the vector is consistent with the smoothing_type parameter. For example, distance_thre = c(NA,round(2/(3^0.5),2),round(4/(3^0.5),2))
#'
#' @return
#' @import tidyverse
#' @export
#'
#' @examples
get_spot_by_celltype_extended=function(
    spatialdf,
    spot_by_celltype,
    smoothing_type = c("loop0","loop1", "loop1+loop2"),
    distance_thre
){

  spatialdf = as.matrix(spatialdf)
  spot_by_celltype = as.matrix(spot_by_celltype)
  smoothing_type = sort(smoothing_type)

  if (identical(smoothing_type,c("loop0","loop1", "loop1+loop2"))) {
    #Make a distance matrix
    distance_matrix = make_distance_matrix(spatialdf = spatialdf,dm_type = setdiff(smoothing_type,"loop0"),distance_thre = distance_thre[-1])

    dst1 = distance_matrix[[smoothing_type[2]]]
    dst2 = distance_matrix[[smoothing_type[3]]]

    spot_by_celltype.dst1=dst1 %*% spot_by_celltype
    colnames(spot_by_celltype.dst1)=paste0(colnames(spot_by_celltype.dst1),paste0("_",smoothing_type[2]))

    spot_by_celltype.dst2=dst2 %*% spot_by_celltype
    colnames(spot_by_celltype.dst2)=paste0(colnames(spot_by_celltype.dst2),paste0("_",smoothing_type[3]))

    if (
      identical(rownames(spot_by_celltype),rownames(spot_by_celltype.dst1)) &
      identical(rownames(spot_by_celltype),rownames(spot_by_celltype.dst2))
    ) {
      spot.ct.df=spot_by_celltype %>% cbind(spot_by_celltype.dst1) %>% cbind(spot_by_celltype.dst2)
      #spot.ct.df
    }
  }

  if (identical(smoothing_type,c("loop0","loop1"))) {
    #Make a distance matrix
    distance_matrix = make_distance_matrix(spatialdf = spatialdf,dm_type = setdiff(smoothing_type,"loop0"),distance_thre = distance_thre[-1])

    dst1 = distance_matrix[[smoothing_type[2]]]
    spot_by_celltype.dst1=dst1 %*% spot_by_celltype
    colnames(spot_by_celltype.dst1)=paste0(colnames(spot_by_celltype.dst1),paste0("_",smoothing_type[2]))

    if (
      identical(rownames(spot_by_celltype),rownames(spot_by_celltype.dst1))
    ) {
      spot.ct.df=spot_by_celltype %>% cbind(spot_by_celltype.dst1)
      #spot.ct.df
    }
  }

  if (identical(smoothing_type,c("loop0"))) {
    spot.ct.df = spot_by_celltype
    #spot.ct.df
  }
  spot.ct.df
}
