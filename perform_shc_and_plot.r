perform_shc_and_plot <- function(sobj, batch = NULL, alpha = 0.25, num_features = 2000, num_PCs = 20, parallel = TRUE, cores = 6) {
  # 检查sobj是否是Seurat对象
  if (!inherits(sobj, "Seurat")) {
    stop("输入的sobj必须是一个Seurat对象")
  }
  
  # 检查group和date是否为字符型
  if (!is.character(group) ||!is.character(date)) {
    stop("group和date必须是字符型")
  }
  
  # 执行scSHC进行显著性检验
  tryCatch({
    shc <- scSHC(GetAssayData(sobj, slot = "data"), batch = NULL, alpha = alpha, num_features = num_features, num_PCs = num_PCs, parallel = parallel, cores = cores)
  }, error = function(e) {
    stop(paste("执行scSHC时出错:", conditionMessage(e)))
  })
  
  # 添加元数据到sobj
  sobj_shc <- AddMetaData(sobj, data.frame(shc[[1]]))
  
  # 绘制UMAP图
  tryCatch({
    p_umap_cluster <- DimPlot(sobj_shc, group.by = "shc..1..", reduction = "umap", pt.size = 0.1) +
      ggtitle("") +
      scale_color_brewer(palette = 'Paired') +
      ggthemes::theme_few() +
      theme(axis.text.x = element_text(size = 15, color = "black"),
            axis.text.y = element_text(size = 15, color = "black"),
            axis.title.x = element_text(size = 15, color = "black"),
            axis.title.y = element_text(size = 15, color = "black"),
            axis.ticks.length.x = unit(0.2, "cm"),
            axis.ticks.length.y = unit(0.2, "cm"),
            axis.ticks.x = element_line(color = "black", size = 1, lineend = 8),
            axis.ticks.y = element_line(color = "black", size = 1, lineend = 8)) +
      theme(plot.title = element_text(hjust = 0.5))
  }, error = function(e) {
    stop(paste("绘制UMAP图时出错:", conditionMessage(e)))
  })
  
  # 保存绘图
  tryCatch({
    ggsave(filename = paste0("dimplot-umap-", "_shc005", ".pdf"), plot = p_umap_cluster, width = 4, height = 4)
  }, error = function(e) {
    warning(paste("保存绘图时出错:", conditionMessage(e)))
  })
  
  # 保存Seurat对象
  tryCatch({
    saveRDS(sobj_shc, paste0("figs/", date, group, "-dimplot-umap-", "_shc005", ".rds"))
  }, error = function(e) {
    warning(paste("保存Seurat对象时出错:", conditionMessage(e)))
  })
  
  return(sobj_shc)
}
