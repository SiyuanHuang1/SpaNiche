#' distribution_2d
#'
#' @param feature_matrix A matrix representing features, which could be topics, interactions, or gene expressions. The matrix should have features as its rows and samples/observations as its columns.
#' @param one_feature A single feature name from the feature_matrix that the user wants to visualize.
#' @param coordinate A data frame or matrix ideally sourced from a Seurat object. It should contain three columns: "row", "col", and "SB".
#' @param thre_to_zero A threshold value used to truncate low expression values.
#' @param color_palette A character vector specifying the color scheme for visualization. The default colors are light grey and blue. If given a single value, it should be one of the following: "viridis", "magma", "cividis", "inferno", "plasma". If given two values, they specify the gradient colors for low and high expression values, respectively.
#' @param visualization_type A character string indicating the type of visualization. The default is "point", but "density" can also be chosen to visualize data density.
#' @param pt.size A numeric value determining the size of the points in the plot.
#' @param axis.title.size A numeric value determining the font size of the axis titles.
#' @param axis.text.size A numeric value determining the font size of the axis text.
#' @param plot.title.size A numeric value determining the font size of the plot's title.
#' @param raster A boolean value indicating whether to rasterize the output plot. This can be useful for plots with a large number of points to improve rendering performance. By default, it's set to FALSE.
#'
#' @return a ggplot2 object
#' @import tidyverse
#' @import Nebulosa
#' @export
#'
#' @examples
distribution_2d=function(
    feature_matrix,
    one_feature,
    coordinate,
    thre_to_zero,
    color_palette = c("lightgrey","blue"),
    visualization_type = "point",
    pt.size = 1,
    axis.title.size = 18,
    axis.text.size = 12,
    plot.title.size = 14,
    raster = F
){


  one.feature.exp=feature_matrix[one_feature,] %>% as.data.frame()
  colnames(one.feature.exp) = "exp"
  one.feature.exp$SB=rownames(one.feature.exp)


  one.feature.exp=one.feature.exp%>%inner_join(coordinate,by = "SB")
  ### 归0
  if (thre_to_zero > 0 & thre_to_zero < 1) {
    thre1 = quantile(one.feature.exp$exp,thre_to_zero)
    one.feature.exp$exp[one.feature.exp$exp < thre1] = 0
  } else {
    stop("range for 'thre_to_zero' can only be (0,1).")
  }


  if(visualization_type == "density"){
    z=Nebulosa:::calculate_density(
      w = one.feature.exp[,1],
      x = one.feature.exp[,3:4],
      method = "wkde",
      adjust = 1, map = TRUE
    )

    colnames(one.feature.exp)[1] = "original_exp"
    one.feature.exp$exp = z
  }


  p=ggplot(data = one.feature.exp,aes(x=col,y=row,color=exp))+
    geom_point(size=pt.size,shape = 16)+
    labs(title = one_feature)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),

      axis.title = element_text(size = axis.title.size),
      axis.text = element_text(color="black",size = axis.text.size),
      axis.ticks.length = unit(0.15,"cm"),

      plot.title = element_text(hjust = 0.5,size = plot.title.size)
    )


  if (length(color_palette) == 1) {
    p=p+scale_color_viridis_c(option = color_palette)
  } else if (length(color_palette) == 2) {
    p=p+scale_color_gradient(low = color_palette[1],high = color_palette[2])
  }

  if(raster){
    ggrastr::rasterise(p, dpi = 300)
  }else{
    p
  }
}
