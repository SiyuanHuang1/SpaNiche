#' spatme_nmf_celltype_and_gene
#'
#' @param vis.seu vis.seu
#' @param spatialdf spatialdf
#' @param enrichmentdf enrichmentdf
#' @param smoothing_type smoothing_type
#' @param distance_thre distance_thre
#' @param ngenes ngenes
#' @param topic_num topic_num
#' @param defined_weight defined_weight
#' @param lambda_v lambda_v
#' @param sigma_v sigma_v
#' @param maxiter maxiter
#' @param st.count st.count
#' @param epsilon epsilon
#'
#' @return
#' @import Seurat
#' @import tidyverse
#' @import IntNMF
#' @import scales
#' @export
#'
#' @examples
spatme_nmf_celltype_and_gene = function(
  ## 基本输入
  vis.seu,
  spatialdf,
  enrichmentdf,
  smoothing_type = c("loop0", "loop1", "loop1+loop2"),
  distance_thre = c(NA,round(2/(3^0.5),2),round(4/(3^0.5),2)),
  ## gene相关参数
  # 可能需要提前选择一些基因
  ngenes = 2000,
  ## 联合分解相关的参数
  topic_num = 15,
  defined_weight = "default", #defined_weight == "default"默认取值；defined_weight = c(0.7,0.3)
  lambda_v =c(0.5, 1, 2),
  sigma_v = c(0.5, 1, 1.5),
  maxiter = 200,
  st.count = 20,
  epsilon = 1e-04
){
  #调整顺序
  if(identical(rownames(spatialdf),rownames(enrichmentdf))){
    print("The spot order is consistent between spatialdf and enrichmentdf.")
  } else {
    enrichmentdf=enrichmentdf[rownames(spatialdf),]
    print("The spot order between spatialdf and enrichmentdf is inconsistent. Adjust the spot order of enrichmentdf.")
  }

  #拓展矩阵
  spot.ct.df = get_spot_by_celltype_extended(
    spatialdf = spatialdf,
    spot_by_celltype = enrichmentdf,
    smoothing_type = smoothing_type,
    distance_thre = distance_thre
  )
  write.csv(spot.ct.df,file = "spot_by_celltype_extended.csv",quote = F,row.names = T)

  ### 第二部分 ###################################################################
  vis.seu <- NormalizeData(vis.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  vis.seu <- FindVariableFeatures(vis.seu, selection.method = "vst", nfeatures = ngenes)
  genedf <- vis.seu@assays$Spatial@data[VariableFeatures(vis.seu),] %>% as.matrix() %>% t()

  #调整顺序
  if(identical(rownames(spatialdf),rownames(genedf))){
    print("The spot order is consistent between spatialdf and genedf.")
  } else {
    genedf=genedf[rownames(spatialdf),]
    print("The spot order between spatialdf and genedf is inconsistent. Adjust the spot order of genedf.")
  }

  #拓展矩阵
  spot.gene.df = get_spot_by_celltype_extended(
    spatialdf = spatialdf,
    spot_by_celltype = genedf,
    smoothing_type = smoothing_type,
    distance_thre = distance_thre
  )

  ### 第三部分 ###################################################################

  #转成matrix
  spot.ct.df = as.matrix(spot.ct.df)
  spot.gene.df = as.matrix(spot.gene.df)
  # 如何在第一次跑IntNMF之前，估计较好的权重
  dat <- list(spot.ct.df,spot.gene.df)
  theta_v = define_initial_weight(dat)
  # The function nmf.mnnals requires the samples to be on rows and variables on columns.
  fit <- nmf.mnnals(dat=dat,k=topic_num,maxiter=50,st.count=10,n.ini=3,ini.nndsvd=TRUE,seed=TRUE,wt = theta_v)
  # 如何确定更严谨的权重
  if (length(defined_weight) == 2 & all(defined_weight > 0)) {
    rho_v = define_good_weight(X = dat,Wzero = fit$W,Hzero = fit$H,theta_v = theta_v,weight_input = defined_weight)
  } else if (length(defined_weight) == 1 & defined_weight == "default") {
    rho_v = theta_v
  } else {
    stop("defined_weight error")
  }
  print("Weights for different matrices have been determined, now starting joint factorization of the two matrices:")

  # 联合分解
  myfit = nmf_modified(
    dat = dat,
    Wzero = fit$W,
    Hzero = fit$H,
    wt = rho_v,
    k=topic_num,
    lambda = lambda_v,
    sigma = sigma_v,
    maxiter = maxiter,
    st.count = st.count,
    spatial.mat = as.matrix(spatialdf[,c("row","col")]),
    epsilon = epsilon
  )
  saveRDS(myfit,file = "nmf_modified_fit.rds")
  print("joint factorization has been completed.")


  output = list(
    #LRintegratedmatrix=LRintegratedmatrix,
    nmf_res=myfit
  )

  print("Output all results.")
  return(output)
}
