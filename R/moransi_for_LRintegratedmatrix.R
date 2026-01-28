#' Calculate Moran's I for spatial ligand-receptor integrated matrix
#'
#' @param LRintegratedmatrix A numeric matrix of integrated LR interaction
#'   scores. Rows correspond to LR interactions and columns correspond to
#'   spatial spots. Row names represent interaction identifiers.
#'   Underscores ('_') in row names will be replaced by hyphens ('-') to ensure
#'   consistency with Seurat's feature naming conventions.
#' @param spatial.seu A Seurat object containing spatial transcriptomics data.
#'   This object must include spatial image information (e.g. Visium) and
#'   spot coordinates. Column names of LRintegratedmatrix must match
#'   the spot names in this object.
#' @param cc_interaction A data frame containing ligand–receptor interaction
#'   annotations.
#'
#' @return A data frame containing Moran's I statistics for spatially variable
#'   ligand–receptor interactions, merged with ligand–receptor annotation
#'   information. The output is sorted by decreasing observed Moran's I.
#' @import tidyverse
#' @import Seurat
#' @export
moransi_for_LRintegratedmatrix=function(
    LRintegratedmatrix,
    spatial.seu,
    cc_interaction
){
  print("calculating moransi for LR integrated matrix:")

  ### Calculate Moran's I using Seurat
  LRintegratedmatrix=LRintegratedmatrix[rowSums(LRintegratedmatrix) > 0,]
  # Replace underscores with hyphens to match Seurat's feature naming convention
  rownames(LRintegratedmatrix)=rownames(LRintegratedmatrix) %>% str_replace_all("_","-")
  # if (!identical(colnames(LRintegratedmatrix), colnames(spatial.seu))) {
  #   stop("Column names or order in LRintegratedmatrix are not identical to spatial.seu!")
  # }
  # Ensure that all spatial spots are present in the LR matrix
  stopifnot(all(colnames(spatial.seu) %in% colnames(LRintegratedmatrix)))
  # Reorder LR matrix columns to match the spatial Seurat object
  LRintegratedmatrix=LRintegratedmatrix[,colnames(spatial.seu)]
  # Create a Seurat object treating LR interactions as features
  # Feature names have already been harmonized to avoid automatic renaming
  pair.seu=CreateSeuratObject(LRintegratedmatrix,assay = "Spatial")
  pair.seu[["Spatial"]]@scale.data = as.matrix(LRintegratedmatrix)
  # Copy spatial image information from the original spatial Seurat object
  pair.seu@images=spatial.seu@images


  # Identify spatially variable LR interactions using Moran's I
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
  # Ensure consistent ligand->receptor direction
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

  ### Merge Moran's I results with LR annotation table
  moransi_output_df = cc_interaction %>% inner_join(moransi_output_df,by = "interaction_name_new")
  print("Finished.")
  print("------------------------------------")
  return(moransi_output_df)
}
