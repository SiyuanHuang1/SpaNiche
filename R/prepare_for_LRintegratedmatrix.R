#' prepare_for_LRintegratedmatrix
#'
#' @param seu_obj ST data is processed as a Seurat object with 'Spatial' as the assay name. Seurat object is normalized by using the 'NormalizeData' function.
#' @param ccc_db This parameter represents the reference for cell-cell interactions, which comes from CellChat.
#'
#' @return
#' @import tidyverse
#' @import Seurat
#' @export
#'
#' @examples
prepare_for_LRintegratedmatrix=function(
    seu_obj,
    ccc_db
){
  wehave_genes=rownames(seu_obj)


  ### 仅基于ccc_db
  cc_interaction=ccc_db$interaction
  cc_complex=ccc_db$complex
  cc_geneInfo=ccc_db$geneInfo
  # 有多少带complex的受/配体
  complex1=setdiff(union(cc_interaction$ligand, cc_interaction$receptor) ,cc_geneInfo$Symbol)


  ### 借助自己的数据对上述参考做过滤
  cc_geneInfo=cc_geneInfo%>%dplyr::filter(Symbol %in% wehave_genes)
  tmp_complex=cc_complex%>%apply(2,function(x){x %in% c(wehave_genes,"")})
  cc_complex=cc_complex[!rowSums(tmp_complex) < 4,]
  # 过滤后有多少候选的complex，以及它的基因矩阵
  complex1=intersect(complex1,rownames(cc_complex))
  cc_complex=cc_complex[complex1,]


  ### 提取基因
  # cc_complex涉及的基因
  gene1=unlist(cc_complex)
  gene1=gene1[gene1 != ""]
  gene1=as.character(gene1)
  # 普通受配体涉及的基因
  gene2=intersect(union(cc_interaction$ligand, cc_interaction$receptor) ,cc_geneInfo$Symbol)
  # 在基于自己数据的前提下，cc_interaction涉及到的所有基因
  gene.used=union(gene1,gene2)


  ### 提取数据
  gene.used=intersect(rownames(seu_obj[["Spatial"]]@data),gene.used) #我的数据包含的基因，调整顺序
  spatial.data=seu_obj[["Spatial"]]@data[gene.used,] %>% as.data.frame() #我的数据缩小范围


  ### 补充complex表达数据
  new.data=matrix(nrow = length(complex1),ncol = dim(seu_obj)[2]) #complex的表达矩阵

  for (i in 1:length(complex1)) {
    one_complex=cc_complex[complex1[i],] %>% as.character()
    one_complex=one_complex[one_complex!=""]

    one_complex_mean=colMeans(spatial.data[one_complex,])
    new.data[i,] = one_complex_mean
  }
  new.data=as.data.frame(new.data)
  colnames(new.data)=colnames(spatial.data)
  rownames(new.data)=complex1

  new.data=rbind(new.data,spatial.data[gene2,])


  ### 过滤原始的互作df
  cc_interaction=cc_interaction[cc_interaction$ligand %in% rownames(new.data) & cc_interaction$receptor %in% rownames(new.data),]


  ### 导出结果
  list(new.data,cc_interaction)
}
