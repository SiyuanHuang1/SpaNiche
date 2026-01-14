#' make_distance_matrix
#'
#' @description
#' The function calculates Euclidean distances between spatial spots and converts them into adjacency matrices based on specified distance thresholds.
#'
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". Row names should  correspond to spot barcodes.
#' @param dm_type This parameter can be one or more values from c('view1', 'view2'), representing the type of adjacency matrix to be returned. For example, you can choose c('view1') or c('view1', 'view2').
#' @param distance_thre This parameter represents distance, defining the central distance from center to it's first/second-order neighbors. Consistent with dm_type, you can define c(c1) or c(c1, c2), and c1/c2 is distance value.
#' @param digits integer indicating the number of decimal places.
#'
#' @return A named list of adjacency matrices corresponding to each dm_type.
#' Each matrix is normalized so that rows sum to 1 (row-normalization matrix).
#'
#' @import tidyverse
#' @importFrom raster pointDistance
#' @export
#'
#' @examples
make_distance_matrix = function(
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
  names(distance_thre) = dm_type

  dst <- pointDistance(spatialdf,lonlat=F)
  dst=round(dst,digits)

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
