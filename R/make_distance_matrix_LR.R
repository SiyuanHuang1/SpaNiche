#' make_distance_matrix_LR
#'
#' @description
#' This function computes pairwise Euclidean distances between spatial spots and constructs distance-based adjacency matrices designed for ligand-receptor interaction information.
#'
#' @param spatialdf A data frame or matrix containing spatial coordinates. It must include two columns named \code{col} and \code{row}. Row names should correspond to spot barcodes.
#' @param dm_type This parameter can be one or more values from c('view1', 'view2'), representing the type of adjacency matrix to be returned.
#' @param distance_thre A numeric vector specifying distance thresholds. If only \code{"view1"} is requested, a single value is sufficient, defining the distance from center to first-order neighbors (view1). If \code{"view1"} + \code{"view2"} is requested, two values are required, defining the distance from center to first-order neighbors (view1) and from center to second-order neighbors (view2).
#' @param digits integer indicating the number of decimal places.
#'
#' @return A named list of row-normalized adjacency matrices.
#' @import tidyverse
#' @importFrom raster pointDistance
make_distance_matrix_LR = function(
    spatialdf,
    dm_type = c("view1", "view2"),
    distance_thre,
    digits = 2
){
  distance_thre = round(distance_thre,digits)

  if (!all(c("col", "row") %in% colnames(spatialdf))) {
    stop("spatialdf must contain columns named 'col' and 'row'")
  }
  spatialdf=spatialdf[,c("col","row")]
  spatialdf=as.matrix(spatialdf)
  if (length(dm_type) != length(distance_thre)) {
    stop("Length of 'distance_thre' must match 'dm_type'")
  }
  dm_type=sort(dm_type)
  distance_thre=sort(distance_thre)
  # names(distance_thre) = dm_type

  dst <- pointDistance(spatialdf,lonlat=F)
  dst=round(dst,digits)

  res_list=as.list(1:length(dm_type))
  names(res_list)=dm_type

  if ("view1" %in% dm_type) {
    thre1 = distance_thre[1]
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
    thre1 = distance_thre[1]
    thre2 = distance_thre[2]
    dst_view2=dst

    dst_view2[dst_view2 == 0] = 1
    dst_view2[dst_view2 > 0 & dst_view2 <= thre1] = 0
    dst_view2[dst_view2 > thre1 & dst_view2 <= thre2] = 1
    dst_view2[dst_view2 > thre2] = 0
    colnames(dst_view2)=rownames(spatialdf)
    rownames(dst_view2)=rownames(spatialdf)
    coef2=rowSums(dst_view2)
    coef2=coef2 ^ (-1)
    dst_view2=dst_view2 %>% apply(2,function(x){ x*coef2})
    res_list$view2 = dst_view2
  }

  res_list
}
