integrate_modalities <- function(rna, m6A) {
  # 共享细胞筛选
  common_cells <- intersect(colnames(rna), colnames(m6A))
  
  # 元数据整合
  rna <- AddMetaData(rna, m6A@meta.data[common_cells, ])
  m6A <- AddMetaData(m6A, rna@meta.data[common_cells, ])
  
  # 创建多组学对象
  multi_obj <- list(
    RNA = CreateSeuratObject(rna@assays$RNA@counts),
    m6A = CreateSeuratObject(m6A@assays$peaks@counts)
  )
  
  return(multi_obj)
}
