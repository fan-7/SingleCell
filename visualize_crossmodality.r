visualize_crossmodality <- function(rna, m6A) {
  # 检查输入对象是否有效
  if (!("Seurat" %in% class(rna)) ||!("Seurat" %in% class(m6A))) {
    stop("输入的对象不是Seurat对象")
  }
  if (!("umap" %in% names(rna@reductions)) ||!("umap" %in% names(m6A@reductions))) {
    stop("输入的Seurat对象中不存在umap降维结果")
  }
  
  # UMAP坐标对齐
  rna$modality <- "RNA"
  m6A$modality <- "m6A"
  seurat <- list(rna, m6A)
  
  nmodalities      = 2
  umap_embeddings <- list()
  
  for(i in 1:length(seurat[1:nmodalities])){
    umap_embeddings[[i]]               <- as.data.frame(seurat[[i]]@reductions[['umap']]@cell.embeddings)
    umap_embeddings[[i]]$UMAP_1        <- umap_embeddings[[i]]$UMAP_1 + (i-1)*40
    umap_embeddings[[i]]$modality      <- unique(seurat[[i]]$modality)
    umap_embeddings[[i]]$cluster       <- seurat[[i]]$seurat_clusters
    umap_embeddings[[i]]$cell_barcode  <- rownames(umap_embeddings[[i]])
  }
  
  umap.embeddings.merge               <- purrr::reduce(umap_embeddings,rbind)
  common.cells                        <- table(umap.embeddings.merge$cell_barcode)
  common.cells                        <- names(common.cells[common.cells==nmodalities])
  
  umap.embeddings.merge               <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]
  
  # 创建连接图
  p <- ggplot(data = umap.embeddings.merge, aes(x = UMAP_1, y = UMAP_2, col = cluster)) +
    geom_point(size = 0.2) +
    scale_color_brewer(palette = 'Pastel1') +
    geom_line(data = umap.embeddings.merge, aes(group = cell_barcode, col = cluster), alpha = 0.1, size = 0.5) +
    theme_classic() + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = "none") +
    labs(title = "Cross-modality Cell Alignment")
  
  return(p)
}
