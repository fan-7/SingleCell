visualize_qc_metrics <- function(sobj, group_name) {
  # 元数据处理
  
  meta_df <- sobj@meta.data %>%
    mutate(UMIs = nCount_RNA,
           Genes = nFeature_RNA,
           UmisPerGene = UMIs/Genes) %>%
    pivot_longer(cols = c(UMIs, Genes), 
                 names_to = "Type", 
                 values_to = "Number")
  
  # 颜色方案
  HMcolors <- c("#D6604D", "#4393C3")
  names(HMcolors) <- c("UMIs","Genes")
  
  # 绘图主题
  my_theme <- theme(
    axis.text = element_text(size=12, color="black"),
    axis.title = element_text(size=14, face="bold"),
    panel.border = element_rect(fill=NA, color="black", size=1),
    panel.grid = element_blank(),
    legend.position = "top"
  )
  
  # 密度图
  density_plot <- ggplot(meta_df, aes(x=Number, fill=Type)) +
    geom_density(alpha=0.7) +
    scale_x_log10() +
    scale_fill_manual(values=HMcolors) +
    labs(title=paste("QC Metrics -", group_name)) +
    my_theme
  # 箱线图
  box_plot <- ggplot(meta_df, aes(x=Type, y=log10(Number), fill=Type)) +
    geom_violin(alpha=0.7) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    scale_fill_manual(values=HMcolors) +
    my_theme
  
  # 组合输出
  gridExtra::grid.arrange(density_plot, box_plot, ncol=2)
}
