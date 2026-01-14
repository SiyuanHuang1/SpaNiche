#' make_distance_matrix
#'
#' @description
#' The function calculates Euclidean distances between spatial spots and converts them into adjacency matrices based on specified distance thresholds.
#'
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". The rownames of this data frame represent the barcode of the spot.
#' @param dm_type This parameter can be one or more values from c('view1', 'view2'), representing the type of adjacency matrix to be returned. For example, you can choose c('view1') or c('view1', 'view2')
#' @param distance_thre This parameter represents distance, defining the central distance from center to it's first/second-order neighbors. consistent with dm_type, you can define c(c1) or c(c1, c2), c1/c2 is distance value.
#'
#' @return A named list of adjacency matrices corresponding to each dm_type.
#' Each matrix is normalized so that rows sum to 1 (row-stochastic matrix).
#' Diagonal elements are 0 (no self-connections).
#'
#' @import tidyverse
#' @importFrom raster pointDistance
#' @export
#'
#' @examples
make_distance_matrix = function(
    spatialdf,
    dm_type = c("view1", "view2"),
    distance_thre
    #distance_thre = c(round(2/(3^0.5),2),round(4/(3^0.5),2))
){

  spatialdf=spatialdf[,c("col","row")]
  spatialdf=as.matrix(spatialdf)
  dm_type=sort(dm_type)
  distance_thre=sort(distance_thre)
  names(distance_thre) = dm_type

  dst <- pointDistance(spatialdf,lonlat=F)
  dst=round(dst,2)

  res_list=as.list(1:length(dm_type))
  names(res_list)=dm_type

  if ("view1" %in% dm_type) {
    thre1 = distance_thre["view1"]
    dst_view1=dst

    dst_view1[dst_view1 <= thre1] = 1
    dst_view1[dst_view1 > thre1] = 0
    colnames(dst_view1)=rownames(spatialdf)
    rownames(dst_view1)=rownames(spatialdf)
    coef1=rowSums(dst_view1)
    coef1=coef1 ^ (-1)
    dst_view1=dst_view1 %>% apply(2,function(x){ x*coef1})
    res_list$view1 = dst_view1
  }

  if ("view2" %in% dm_type) {
    thre2 = distance_thre["view2"]
    dst_view2=dst

    dst_view2[dst_view2 <= thre2] = 1
    dst_view2[dst_view2 > thre2] = 0
    colnames(dst_view2)=rownames(spatialdf)
    rownames(dst_view2)=rownames(spatialdf)
    coef2=rowSums(dst_view2)
    coef2=coef2 ^ (-1)
    dst_view2=dst_view2 %>% apply(2,function(x){ x*coef2})
    res_list[["view2"]] = dst_view2
  }

  res_list
}
