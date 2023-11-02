#' spatme_main
#'
#' @param vis.seu ST data is processed as a Seurat object with 'Spatial' as the assay name. Seurat object is normalized by using the 'NormalizeData' function.
#' @param spatialdf a data frame or matrix that contains two columns, "col" and "row". The row names represent the barcode of the spot.
#' @param enrichmentdf a data frame or matrix representing the abundance of cell types, with spots as rows and cell types as columns. Its row names should be consistent with the row names of 'spatialdf'.
#' @param smoothing_type This parameter has three options: c('loop0'), c('loop0', 'loop1'), and c('loop0', 'loop1', 'loop1+loop2'). 'loop0', 'loop1', and 'loop1+loop2' indicate the range within which data smoothing is performed, denoting: no smoothing, smoothing using one adjacent spot, and smoothing using two surrounding spots, respectively. As it involves the merging of matrices, the returned matrix's dimension varies depending on the chosen option. The dimensions corresponding to the three options are 1x, 2x, and 3x of the original dimensions, respectively.
#' @param distance_thre This parameter represents distance, defining the central distance from loop0 to loop1, and from loop0 to loop2. The length of the vector is consistent with the smoothing_type parameter. For example, distance_thre = c(NA,round(2/(3^0.5),2),round(4/(3^0.5),2))
#' @param LRDB This parameter represents the reference for cell-cell interactions, which comes from CellChat.
#' @param LR_dm_type This parameter can be one or more values from c('loop1', 'loop2'), representing the type of adjacency matrix to be returned.
#' @param LR_distance_thre This parameter represents distance, defining the central distance from loop0 to loop1(, and from loop0 to loop2).
#' @param var_thre Variance threshold
#' @param topic_num An integer specifying the number of components or topics to be extracted. It determines the number of columns in the W matrix and the number of rows in each H matrix.
#' @param defined_weight defined_weight
#' @param lambda_v lambda
#' @param sigma_v Gaussian kernel parameter. Used in the computation of the Laplacian matrix, influencing how spatial information is incorporated into the factorization.
#' @param maxiter The maximum number of iterations allowed for the factorization process. This acts as a stopping criterion to prevent the algorithm from running indefinitely.
#' @param st.count Convergence counter. If the change in the reconstruction error remains below epsilon for st.count consecutive iterations, the algorithm stops, assuming it has converged.
#' @param epsilon A threshold for the change in reconstruction error. If the relative change in error is less than epsilon for st.count consecutive iterations, convergence is assumed.
#' @param topic_lr_quantile_thre A threshold for quantile selection. Only data above this quantile will be considered for visualization. Default is 0.95.
#' @param spot_topic.quantile.thre spot_topic.quantile.thre
#' @param S_intra_thre S_intra_thre
#' @param term_topn term_topn
#' @param num_core Number of core to be used during parallel computation. only be used when running scSeqComm analyze
#'
#' @return
#' @import Seurat
#' @import tidyverse
#' @import IntNMF
#' @import scales
#' @import scSeqComm
#' @export
#'
#' @examples
spatme_main = function(
    ## 基本输入
  vis.seu,
  spatialdf,
  enrichmentdf,
  smoothing_type = c("loop0", "loop1", "loop1+loop2"),
  distance_thre = c(NA,round(2/(3^0.5),2),round(4/(3^0.5),2)),
  ## LR相关参数
  LRDB = spatme::CellChatDB.human,
  LR_dm_type = c("loop1", "loop2"), #暂时不能更改
  LR_distance_thre = c(round(2/(3^0.5),2),round(4/(3^0.5),2)),
  var_thre = c(0.2,0.99),
  ## 联合分解相关的参数
  topic_num = 15,
  defined_weight = "default", #defined_weight == "default"默认取值；defined_weight = c(0.7,0.3)
  lambda_v =c(0.5, 1, 2),
  sigma_v = c(0.5, 1, 1.5),
  maxiter = 200,
  st.count = 20,
  epsilon = 1e-04,
  ## 下游的细胞响应
  topic_lr_quantile_thre = 0.95,
  spot_topic.quantile.thre = 0.2,
  S_intra_thre = 0.8,
  term_topn = 20,
  num_core = 1
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
  step1_res=prepare_for_LRintegratedmatrix(vis.seu,LRDB)
  distance_matrix = make_distance_matrix_LR(spatialdf = spatialdf,dm_type = LR_dm_type,distance_thre = LR_distance_thre)
  LRelementmatrix = step1_res[[1]]
  cc_interaction = step1_res[[2]]

  # 单个spot中，R的表达量
  rmat = LRelementmatrix %>% t() %>% .[,cc_interaction$receptor]
  # 单个spot中，L的表达量
  lmat = LRelementmatrix %>% t() %>% .[,cc_interaction$ligand]
  # 后续为了扩大适用性，下面要改，不一定是m1到m4
  m1=generate_LRintegratedmatrix(neighbor_matrix = distance_matrix$loop1,neighbor_role = rmat,center_role = lmat,center_is_ligand = T)
  m2=generate_LRintegratedmatrix(neighbor_matrix = distance_matrix$loop1,neighbor_role = lmat,center_role = rmat,center_is_ligand = F)
  colnames(m1)=paste0(colnames(m1)," (loop0+loop1)")
  colnames(m2)=paste0(colnames(m2)," (loop0+loop1)")

  m3=generate_LRintegratedmatrix(neighbor_matrix = distance_matrix$loop2,neighbor_role = rmat,center_role = lmat,center_is_ligand = T)
  m4=generate_LRintegratedmatrix(neighbor_matrix = distance_matrix$loop2,neighbor_role = lmat,center_role = rmat,center_is_ligand = F)
  colnames(m3)=paste0(colnames(m3)," (loop0+loop2)")
  colnames(m4)=paste0(colnames(m4)," (loop0+loop2)")

  LRintegratedmatrix=m1 %>% cbind(m2) %>% cbind(m3) %>% cbind(m4)
  LRintegratedmatrix[is.na(LRintegratedmatrix)] = 0
  saveRDS(LRintegratedmatrix,file = "LRintegratedmatrix.rds")
  rm(list = paste0("m",1:4))

  ### 计算莫兰指数，然后添加interaction信息 ###
  moransi_output_df = moransi_for_LRintegratedmatrix(t(LRintegratedmatrix),vis.seu,cc_interaction)
  saveRDS(moransi_output_df,file = "LRintegratedmatrix_with_moransi.rds")

  print("Moran's index calculation completed")

  ### 第三部分 ###################################################################
  LRintegratedmatrix = LRintegratedmatrix[,colSums(LRintegratedmatrix) > 0]
  colnames(LRintegratedmatrix)=colnames(LRintegratedmatrix) %>% str_replace_all("_","-")
  LRintegratedmatrix=LRintegratedmatrix[spot.ct.df %>% rownames(),]

  # 方差情况
  var.df1 = LRintegratedmatrix %>% apply(2, function(x){var(x)}) %>% as.data.frame()
  colnames(var.df1) = "var_value"
  var.df1$interaction = rownames(var.df1)
  var.df1 = var.df1 %>% arrange(var_value)
  var.df1$index = 1:length(var.df1$var_value)
  var_thre_1 = quantile(var.df1$var_value,var_thre[1])
  var_thre_2 = quantile(var.df1$var_value,var_thre[2])
  var.df1 = var.df1 %>% filter(var_value > var_thre_1 & var_value < var_thre_2)
  used_interaction = var.df1$interaction
  used_interaction = intersect(colnames(LRintegratedmatrix),used_interaction)
  LRintegratedmatrix=LRintegratedmatrix[,used_interaction]

  print("Filtering the LR integrated matrix based on variance and selecting some LR pairs.")

  #转成matrix
  spot.ct.df = as.matrix(spot.ct.df)
  LRintegratedmatrix = as.matrix(LRintegratedmatrix)
  # 如何在第一次跑IntNMF之前，估计较好的权重
  dat <- list(spot.ct.df,LRintegratedmatrix)
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

  ### 第五部分 ###################################################################
  # 联系下游细胞响应
  spot_topic = as.data.frame(myfit$W)
  colnames(spot_topic) = paste0("topic",1:topic_num)
  allSB = rownames(spot_topic)

  gene_expr_matrix = vis.seu@assays$Spatial@data
  gene_expr_matrix = as.matrix(gene_expr_matrix)
  gene_expr_matrix = gene_expr_matrix[rowSums(gene_expr_matrix) > 0,]

  tmpplot = plot_topic_lr(
    topic_lr = t(myfit$H$H2),
    topic_lr_quantile_thre = topic_lr_quantile_thre
  )
  topic_lr_small_file = dir(getwd(),"topic_lr_small_.*csv")
  topic_lr_small = read.csv(topic_lr_small_file)

  # 配色
  color_cluster = hue_pal()(topic_num)
  names(color_cluster) = paste0("topic",1:topic_num)

  # 需要导出
  plot.list=list()
  all.go = data.frame()

  print("Exploring pathway changes due to cell interactions:")

  for (ti in colnames(spot_topic)) {
    tmpst = spot_topic[,ti]
    tmpquan = quantile(tmpst,c(spot_topic.quantile.thre,1 - spot_topic.quantile.thre))
    low.index = which(tmpst <= tmpquan[1])
    high.index = which(tmpst > tmpquan[2])
    lowSB = allSB[low.index]
    highSB = allSB[high.index]


    ########## Load scRNA-seq data
    tmp_gene_expr_matrix = gene_expr_matrix[,c(highSB,lowSB)]
    tmp_anno = data.frame(
      "Cell_ID" = c(highSB,lowSB),
      "Cluster_ID" = c(
        rep(paste(ti,"high",sep = "_"),length(highSB)),
        rep(paste(ti,"low",sep = "_"),length(lowSB))
      )
    )
    ########## Ligand-receptor pairs
    tmp_topic_lr = topic_lr_small %>% filter(topic == ti)
    tmp_topic_lr$lr = str_remove(tmp_topic_lr$LRinteraction," .*$")
    tmp_topic_lr$r = ifelse(str_detect(tmp_topic_lr$lr,">"),
                            str_remove(tmp_topic_lr$lr,"^.*>"),
                            str_remove(tmp_topic_lr$lr,"<.*$"))
    tmp_topic_lr$l = ifelse(str_detect(tmp_topic_lr$lr,">"),
                            str_remove(tmp_topic_lr$lr,"->.*$"),
                            str_remove(tmp_topic_lr$lr,"^.*<-"))
    tmp_topic_lr=tmp_topic_lr[,c("l","r")] %>% unique()
    colnames(tmp_topic_lr) = c("ligand","receptor")

    tmp_topic_lr$ligand = str_replace(tmp_topic_lr$ligand,"^HLA-","xxx") %>% str_replace("-",",") %>% str_replace("xxx","HLA-")
    tmp_topic_lr$receptor = str_replace(tmp_topic_lr$receptor,"-",",")
    LR_db = tmp_topic_lr
    ########## Transcriptional regulatory networks
    TF_TG_db <- scSeqComm::TF_TG_TRRUSTv2_HTRIdb_RegNetwork_High
    ########## Receptor-Transcription factor a-priori association
    TF_PPR <- scSeqComm::TF_PPR_KEGG_human


    ########## Identify and quantify intercellular and intracellular signaling
    scSeqComm_res <- scSeqComm_analyze(gene_expr = tmp_gene_expr_matrix,
                                       cell_metadata = tmp_anno,
                                       inter_signaling = F,
                                       LR_pairs_DB = LR_db,
                                       TF_reg_DB = TF_TG_db,
                                       R_TF_association = TF_PPR,
                                       N_cores = num_core,
                                       DEmethod = "wilcoxon")


    ### （受配体对形成后，）后续的细胞响应
    topic_high_comm = dplyr::filter(scSeqComm_res$comm_results,cluster == paste(ti,"high",sep = "_") & S_intra >= S_intra_thre)
    if (dim(topic_high_comm)[1] < 2) {next} #20231018

    # GO analysis of topic_high communication
    geneUniverse <- unique(unlist(scSeqComm_res$TF_reg_DB_scrnaseq))
    cell_functional_response <- scSeqComm_GO_analysis(
      results_signaling = topic_high_comm,
      geneUniverse = geneUniverse,
      method = "general")

    # 画图
    cell_functional_response$pval = as.numeric(cell_functional_response$pval)
    cell_functional_response = cell_functional_response %>% arrange(pval)
    cell_functional_response$pval_log10_neg = -log10(cell_functional_response$pval)
    cell_functional_response$cluster = ti
    all.go = rbind(all.go,cell_functional_response %>% filter(pval < 0.01))
    tmpres = cell_functional_response %>% slice_head(n = term_topn)

    tmpres=tmpres%>%arrange(pval_log10_neg)
    tmpres$Term=factor(tmpres$Term,levels = tmpres$Term)

    tmpbar = tmpres %>% ggplot(aes(x=Term,y=pval_log10_neg))+
      geom_hline(yintercept = -log10(0.01),color = "black",alpha=0.7)+
      geom_bar(stat="identity",alpha=0.8,aes(fill=cluster))+
      geom_text(mapping = aes(x=Term,y=0,label=Term),hjust=0)+
      scale_x_discrete("")+
      scale_y_continuous("-log10(p value)",expand = c(0.02,0))+
      scale_fill_manual(values = color_cluster)+
      coord_flip()+
      labs(title = ti)+
      theme_bw()+
      theme(
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 20)
      )

    index=which(ti == colnames(spot_topic))
    plot.list[[get("index")]]=tmpbar
  }
  print("Completed.")

  output = list(
    LRintegratedmatrix=LRintegratedmatrix,
    nmf_res=myfit,
    go_res=all.go,
    go_plot=plot.list
  )

  print("Output all results.")
  return(output)
}
