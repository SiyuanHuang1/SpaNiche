estimate_ecotype_ratio = function(
    sample_gene,
    sample_pair,
    ecotype_gene,
    ecotype_pair,
    dataset.w = c(4,1)
  ){
  if (
    identical(rownames(sample_pair),rownames(sample_gene)) &
    identical(rownames(ecotype_pair),rownames(ecotype_gene)) &
    identical(colnames(sample_pair),colnames(ecotype_pair)) &
    identical(colnames(sample_gene),colnames(ecotype_gene))
  ) {
    print("[dimensions] and [row/column names] of 4 matrices meet the requirements.")
  } else {
    stop("[dimensions] or [row/column names] of 4 matrices do not meet the requirements.")
  }
  
  bulksamples = rownames(sample_gene)
  loss.shouldbe = dataset.w[1] / dataset.w[2]
  loss.shouldbe = loss.shouldbe * c(0.9,1.1) #loss之比应该落在什么范围
  
  sample_ecotype = matrix(0,nrow = nrow(sample_gene),ncol = nrow(ecotype_gene))
  rownames(sample_ecotype) = rownames(sample_gene)
  colnames(sample_ecotype) = rownames(ecotype_gene)
  
  for (bi in bulksamples) {
    # 不额外设置两套数据的权重
    y1 = sample_gene[bi,]
    y2 = sample_pair[bi,]
    y = as.matrix(c(y1,y2))
    
    x1 = t(ecotype_gene)
    x2 = t(ecotype_pair)
    x = rbind(x1,x2)
    
    tmpa = t(x) %*% x
    tmpb = solve(tmpa)
    tmpc = tmpb %*% t(x)
    tmp.res = tmpc %*% y
    
    loss1 = sum((x1 %*% tmp.res - y1) ^ 2)
    loss2 = sum((x2 %*% tmp.res - y2) ^ 2)
    loss.wehave = loss1 / loss2
    
    # 此时的loss占比是否合适？需要判断一下
    # 如果不合适，继续“打磨”
    counter = 0
    while ((!((loss.wehave >= loss.shouldbe[1]) & (loss.wehave <= loss.shouldbe[2]))) & (counter <= 20)) {
      counter = counter + 1
      
      # 调整权重
      # k1 = dataset.w[1] / loss1
      # k2 = dataset.w[2] / loss2
      k1 = sqrt(dataset.w[1] / loss1)
      k2 = sqrt(dataset.w[2] / loss2)
      y1 = y1 * k1
      y2 = y2 * k2
      x1 = x1 * k1
      x2 = x2 * k2
      x = rbind(x1,x2)
      y = as.matrix(c(y1,y2))
      
      tmp.res = solve(t(x) %*% x) %*% t(x) %*% y
      loss1 = sum((x1 %*% tmp.res - y1) ^ 2)
      loss2 = sum((x2 %*% tmp.res - y2) ^ 2)
      loss.wehave = loss1 / loss2
    }
    
    tmp.res[tmp.res < 0] = 0
    tmp.res[,1] = tmp.res[,1] / sum(tmp.res[,1])
    
    sample_ecotype[bi,] = t(tmp.res)
  }
  sample_ecotype = as.data.frame(sample_ecotype)
  return(sample_ecotype)
}