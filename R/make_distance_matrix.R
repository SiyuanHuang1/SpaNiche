#' make_distance_matrix
#'
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". The row names represent the barcode of the spot.
#' @param dm_type This parameter can be one or more values from c('loop1', 'loop1+loop2'), representing the type of adjacency matrix to be returned.
#' @param distance_thre This parameter represents distance, defining the central distance from loop0 to loop1, and/or from loop0 to loop2.
#'
#' @return
#' @import tidyverse
#' @importFrom raster pointDistance
#' @export
#'
#' @examples
make_distance_matrix = function(
    spatialdf,
    dm_type = c("loop1", "loop1+loop2"),
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

  if ("loop1" %in% dm_type) {
    thre1 = distance_thre["loop1"]
    dst_loop1=dst

    dst_loop1[dst_loop1 <= thre1] = 1
    dst_loop1[dst_loop1 > thre1] = 0
    colnames(dst_loop1)=rownames(spatialdf)
    rownames(dst_loop1)=rownames(spatialdf)
    coef1=rowSums(dst_loop1)
    coef1=coef1 ^ (-1)
    dst_loop1=dst_loop1 %>% apply(2,function(x){ x*coef1})
    res_list$loop1 = dst_loop1
  }

  if ("loop1+loop2" %in% dm_type) {
    thre2 = distance_thre["loop1+loop2"]
    dst_loop12=dst

    dst_loop12[dst_loop12 <= thre2] = 1
    dst_loop12[dst_loop12 > thre2] = 0
    colnames(dst_loop12)=rownames(spatialdf)
    rownames(dst_loop12)=rownames(spatialdf)
    coef12=rowSums(dst_loop12)
    coef12=coef12 ^ (-1)
    dst_loop12=dst_loop12 %>% apply(2,function(x){ x*coef12})
    res_list[["loop1+loop2"]] = dst_loop12
  }

  res_list
}
