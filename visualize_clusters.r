visualize_clusters <- function(sobj, group_by) {
  # 颜色方案
  cluster_colors <- scales::hue_pal()(length(unique(sobj[[group_by]])))
  
  # 可视化函数
  plot_dim <- function(reduction) {
    DimPlot(sobj, reduction=reduction, group.by=group_by, pt.size=0.5) +
      scale_color_manual(values=cluster_colors) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  }
  
  # 组合可视化
  plot_umap <- plot_dim("umap")
  plot_tsne <- plot_dim("tsne")
  
  gridExtra::grid.arrange(plot_umap, plot_tsne, ncol=2)
}
