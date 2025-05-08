cytotrace_analysis <- function(sobj) {
  # 从Seurat对象中获取counts数据
  counts <- GetAssayData(object = sobj, slot = "counts")
  
  # 找出counts矩阵中大于0的元素，得到逻辑向量
  nonzero <- counts > 0
  
  # 对逻辑向量按行求和，筛选出在至少10个细胞中表达的基因
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  # 对逻辑向量按列求和，筛选出在至少10个基因中表达的细胞
  keep_cells <- Matrix::colSums(nonzero) >= 10
  
  # 过滤数据，只保留符合条件的基因和细胞
  filtered_counts <- data.frame(counts[keep_genes, keep_cells])
  
  # 修改过滤后数据的行名，这里使用make.names函数确保行名唯一
  rownames(filtered_counts) <- make.names(sapply(strsplit(rownames(filtered_counts), "\\."), "[[", 2), unique = T)
  
  # 运行CytoTRACE分析
  results <- CytoTRACE(filtered_counts)
  
  # 保存分析结果到文件
  saveRDS(results, file = "20240719_result_cytotrace_rna_zlz_E85E105.rds")
  
  # 返回分析结果
  return(results)
}
