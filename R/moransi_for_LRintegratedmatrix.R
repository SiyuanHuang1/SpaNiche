#' moransi_for_LRintegratedmatrix
#'
#' @param LRintegratedmatrix integrated LR matrix, where the row names represent features and the column names represent spots. The underscore ('_') in the row names will be replaced by a hyphen ('-').
#' @param spatial.seu spatial.seu
#' @param cc_interaction cc_interaction
#'
#' @return
#' @import tidyverse
#' @import Seurat
#' @export
#'
#' @examples
moransi_for_LRintegratedmatrix=function(
    LRintegratedmatrix,
    spatial.seu,
    cc_interaction
){
  print("calculating moransi for LR integrated matrix:")

  ### 求莫兰指数
  #借助Seurat
  LRintegratedmatrix=LRintegratedmatrix[rowSums(LRintegratedmatrix) > 0,]
  rownames(LRintegratedmatrix)=rownames(LRintegratedmatrix) %>% str_replace_all("_","-")
  # if (!identical(colnames(LRintegratedmatrix), colnames(spatial.seu))) {
  #   stop("Column names or order in LRintegratedmatrix are not identical to spatial.seu!")
  # }
  LRintegratedmatrix=LRintegratedmatrix[,colnames(spatial.seu)]
  pair.seu=CreateSeuratObject(LRintegratedmatrix,assay = "Spatial") #这一步基因名会变：('_'), replacing with dashes ('-')；提前自己改过来
  pair.seu[["Spatial"]]@scale.data = as.matrix(LRintegratedmatrix)

  pair.seu@images=spatial.seu@images

  pair_moransi <- FindSpatiallyVariableFeatures(
    pair.seu, assay = "Spatial", slot = "scale.data",
    features = rownames(pair.seu),
    selection.method = "moransi"
  )
  moransi_output_df <- pair_moransi@assays$Spatial@meta.features
  moransi_output_df = moransi_output_df %>% arrange(desc(MoransI_observed))
  moransi_output_df$interaction_name_detailed=rownames(moransi_output_df)

  ###
  cc_interaction$interaction_name_new = paste(
    cc_interaction$ligand,cc_interaction$receptor,sep = "->"
    ) %>% str_replace_all("_","-")

  moransi_output_df$interaction_name_new=str_replace(moransi_output_df$interaction_name_detailed," \\(.*","")
  moransi_output_df$interaction_name_new=moransi_output_df$interaction_name_new %>%
    sapply(
      function(x){
        if(str_detect(x,"<-")){
          tmp1=str_split(x,"<-")[[1]][1]
          tmp2=str_split(x,"<-")[[1]][2]
          tmpnew=paste(tmp2,tmp1,sep = "->")
        }else{
          tmpnew=x
        }
        as.character(tmpnew)
      }
    )

  ###
  moransi_output_df = cc_interaction %>% inner_join(moransi_output_df,by = "interaction_name_new")
  print("Finished.")
  print("------------------------------------")
  return(moransi_output_df)
}
