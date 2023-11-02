#' DEinteractions_twoclusters
#'
#' @param spot_topic A data frame where columns represent different topics and row names represent different spots.
#' @param overall.mat An interaction matrix where row names represent interactions and column names represent spots.
#' @param thre.var Variance threshold used to filter interactions.
#' @param thre.topic Topic threshold used to distinguish between the higher and lower clusters.
#' @param thre.p.adj Adjusted p-value threshold. Interactions with an adjusted p-value below this threshold are considered significant.
#' @param thre.log2FC Log2 Fold Change threshold.
#' @param thre.pct.1 Threshold for the proportion of spots in the "high" cluster that have non-zero interaction values.
#' @param thre.pct.2 Threshold for the proportion of spots in the "low" cluster that have non-zero interaction values.
#'
#' @return
#' @import tidyverse
#' @export
#'
#' @examples
DEinteractions_twoclusters = function(
    spot_topic,
    overall.mat,
    thre.var = 0.2,
    thre.topic = 0.9,
    thre.p.adj = 0.01,
    thre.log2FC = 1.5,
    thre.pct.1 = 0.2,
    thre.pct.2 = 0.2
){

  ### 是否顺序一致
  topics = colnames(spot_topic)
  spot_topic$SB=rownames(spot_topic)
  if (identical(rownames(spot_topic),colnames(overall.mat))) {
  } else {
    sameSB=intersect(rownames(spot_topic),colnames(overall.mat))
    overall.mat = overall.mat[,sameSB]
    spot_topic = spot_topic[sameSB,]
  }


  ### 是否可以根据方差先筛一遍
  # 可以。方差很小的interaction不太可能是差异的（不会被找出来）
  interaction.var=overall.mat %>% apply(1, var) %>% as.data.frame()
  colnames(interaction.var)="var"
  interaction.var$interaction=rownames(interaction.var)
  interaction.var=interaction.var%>%arrange(var)
  thre1=quantile(interaction.var$var,thre.var)
  interaction.var=interaction.var%>%dplyr::filter(var >= thre1)

  someinteraction=intersect(rownames(overall.mat),interaction.var$interaction)
  overall.mat=overall.mat[someinteraction,]
  interaction.len=dim(overall.mat)[1]

  ### 开始找差异pair
  all.topic.df=data.frame()

  for (ti in topics) {
    tmpdf=spot_topic[,c("SB",ti)]
    colnames(tmpdf)[2] = "onetopic"
    thre2=quantile(tmpdf$onetopic,thre.topic)

    lowindex=which(tmpdf$onetopic < thre2)
    highindex=which(tmpdf$onetopic >= thre2)

    ### 计算p值和log2FC
    deg.res <- apply(overall.mat, 1, function(x){

      tmpv1=wilcox.test(
        x[highindex],
        x[lowindex]
      )$p.value

      tmpv2=log2(
        (mean(x[highindex])+0.001) /
          (mean(x[lowindex])+0.001)
      )

      tmpv3=sum(x[highindex] > 0) / length(highindex)
      tmpv4=sum(x[lowindex] > 0) / length(lowindex)

      c(tmpv1,tmpv2,tmpv3,tmpv4)
    })
    deg.res=as.data.frame(t(deg.res))
    colnames(deg.res)=c("p.value","log2FC","pct.1","pct.2")
    deg.res$interaction=rownames(deg.res)
    deg.res=deg.res%>%arrange(p.value)

    ### BH
    fdr=rep(0,interaction.len)
    for (i in interaction.len:1) {
      if (i==interaction.len) {
        tmpfdr=deg.res$p.value[i]
      }else{
        tmpfdr=min(tmpfdr,deg.res$p.value[i] * (interaction.len / i))
      }
      fdr[i] = tmpfdr
    }

    deg.res$p.adj=fdr
    deg.res$topic=ti
    all.topic.df=all.topic.df%>%rbind(deg.res)
    print(paste0(ti," is ok!"))
  }

  ### 过滤
  all.topic.df.f=all.topic.df%>%dplyr::filter(
    p.adj < thre.p.adj &
      log2FC >= thre.log2FC &
      pct.1 >= thre.pct.1 &
      pct.2 < thre.pct.2
  )

  list(all.topic.df,all.topic.df.f)
}
