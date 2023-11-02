#' plot_NMF_components
#'
#' @param components A dataframe containing the NMF component data.
#' @param one_topic A string that specifies which topic to visualize. Default is "topic1".
#' @param p1.axis.text.x.bottom.size p1.axis.text.x.bottom.size
#' @param p2.axis.text.x.bottom.size p2.axis.text.x.bottom.size
#' @param p1.strip.text.x.size p1.strip.text.x.size
#' @param p2.strip.text.x.size p2.strip.text.x.size
#' @param color_celltype Either NULL or a vector of colors. If provided, it specifies the colors to use for the different cell types in the plot. If NULL, default ggplot2 colors are used.
#' @param color_loop Either NULL or a vector of colors. If provided, it specifies the colors to use for the different loops in the plot. If NULL, default ggplot2 colors are used.
#'
#' @return A combined ggplot2 plot that visualizes the values of the specified topic grouped by celltype and loop.
#' @import tidyverse
#' @import patchwork
#' @import reshape2
#' @import scales
#' @import RColorBrewer
#' @export
#'
#' @examples
plot_NMF_components = function(
    components,
    one_topic = "topic1",
    p1.axis.text.x.bottom.size = 8,
    p2.axis.text.x.bottom.size = 9,
    p1.strip.text.x.size = 15,
    p2.strip.text.x.size = 9,
    color_celltype = NULL,
    color_loop = NULL
){

  components$topic=rownames(components)
  components=components%>%reshape2::melt(id.vars=c("topic"))
  colnames(components)[2]="celltype_loop"
  components$celltype_loop = as.character(components$celltype_loop)

  components$loop = str_extract(components$celltype_loop,"loop[0-9].*$")
  components$loop[is.na(components$loop)] = "loop0"
  components$celltype=str_replace(components$celltype_loop,"_loop[0-9].*$","")

  each=components%>%dplyr::filter(topic == one_topic)
  ###
  p1=each%>%ggplot(aes(x=celltype,y=value,fill=celltype))+
    geom_bar(stat = "identity")+
    scale_y_continuous(expand = c(0.02,0))+
    facet_grid(.~loop)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,color = "black",size = p1.axis.text.x.bottom.size),
      axis.text.y.left = element_text(color = "black",size = 10),
      axis.title = element_text(size = 15),
      strip.text.x = element_text(size = p1.strip.text.x.size),

      legend.position = "top"
    )

  ###
  p2=each%>%ggplot(aes(x=loop,y=value,fill=loop))+
    geom_bar(stat = "identity")+
    scale_y_continuous(expand = c(0.02,0))+
    facet_grid(.~celltype)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,color = "black",size = p2.axis.text.x.bottom.size),
      axis.text.y.left = element_text(color = "black",size = 10),
      axis.title = element_text(size = 15),
      strip.text.x = element_text(size = p2.strip.text.x.size),

      legend.position = "top"
    )

  ###
  if (is.null(color_celltype)) {
  } else {
    p1 = p1+scale_fill_manual(values = color_celltype)
  }
  if (is.null(color_loop)) {
  } else {
    p2 = p2+scale_fill_manual(values = color_loop)
  }

  ###
  #p1+p2+plot_layout(widths = c(1, relative_width))
  p1 / p2
}
