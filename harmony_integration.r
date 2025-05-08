harmony_integration <- function(sobj, nclust = 50, max_iter = 10, early_stop = TRUE, dims = 1:20, resolution_values = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)) {
  # 执行Harmony整合
  harm <- sobj %>%
    RunHarmony("sample", plot_convergence = TRUE,
               nclust = nclust, max_iter = max_iter, early_stop = early_stop, reduction = "pca", assay.use = "data", reduction.save = "harmony")
  
  # 运行UMAP降维
  harmonized_seurat <- RunUMAP(harm, reduction = "harmony", assay = "data", dims = dims)
  
  # 寻找邻居
  harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
  
  # 执行聚类
  harmonized_seurat <- FindClusters(harmonized_seurat, resolution = resolution_values)
  
  # 绘制DimPlot
  p <- DimPlot(harmonized_seurat, reduction = "harmony")
  
  # 返回整合后的Seurat对象和绘图结果
  return(list(seurat_object = harmonized_seurat, plot = p))
}

# 示例使用
# 假设sobj是一个Seurat对象
# sobj <- readRDS("your_seurat_object.rds")
# result <- harmony_integration(sobj)
# 可以通过以下方式访问整合后的Seurat对象和绘图
# integrated_sobj <- result$seurat_object
# plot_result <- result$plot
# print(plot_result)
