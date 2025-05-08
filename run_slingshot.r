library(slingshot)
library(viridis)

# 定义轨迹推断可视化的函数
run_slingshot_analysis <- function(sobj) {
  sobj$cell_type <- paste0("CT", sobj$seurat_clusters)
  # slingshot推断轨迹 基于UMAP降维
  slingshot.result <- slingshot(data = sobj@reductions[['umap']]@cell.embeddings[, 1:2], 
                                clusterLabels = sobj$cell_type)
  
  # 伪时间分析
  pt <- slingPseudotime(slingshot.result)
  # slingshot结果UMAP降维展示 颜色按伪时间渐变
  cell_colors <- viridis(100, option = "B")[cut(pt[, 1], breaks = 100)]
  plot(slingshot.result$reducedDim,col = cell_colors,pch=16,cex=0.5,axes=FALSE,main='clusters')
  lines(SlingshotDataSet(slingshot.result), lwd = 2,  col = '#454042')
  return(list(slingshot_result = slingshot.result, pseudotime = pt))
}
