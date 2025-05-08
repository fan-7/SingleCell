nmf_decomposition <- function(zlz_cm, ip_cm, rank = 3) {
  # 找到两个矩阵行名的交集，用于后续数据整合
  genes <- intersect(rownames(zlz_cm), rownames(ip_cm))
  
  # 行拼接两个矩阵前，修改矩阵行名
  rownames(zlz_cm) <- paste0("input.", rownames(zlz_cm))
  rownames(ip_cm) <- paste0("ip.", rownames(ip_cm))
  
  # 行拼接两个矩阵
  M <- rbind(zlz_cm[, colnames(ip_cm)], ip_cm)
  
  # 找出矩阵中大于0的元素，得到逻辑向量
  nonzero <- M > 0
  
  # 对逻辑向量按行求和，筛选出在至少10个细胞中表达的基因
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  keep_cells <- Matrix::colSums(nonzero) >= 10
  
  # 过滤数据，只保留符合条件的基因和细胞
  filtered_counts <- data.frame(M[keep_genes, keep_cells])
  
  # 进行NMF分解
  nmf_result <- nmf(filtered_counts, rank, method = "snmf/r", seed = 'nndsvd')
  
  # 提取分解后的矩阵W和H
  W <- coef(nmf_result)
  H <- fit(nmf_result)
  
  # 返回分解后的矩阵W和H
  return(list(W = W, H = H))
}

# 示例使用
zlz_cm <- readRDS("matchIP_cm.rds")
ip_cm <- readRDS("matchzlz_cm.rds")

# 调用NMF函数进行分解
result <- nmf_decomposition(zlz_cm, ip_cm, rank = 3)

# 查看分解后的矩阵W和H 用于下游将为聚类分析
print(result$W)
print(result$H)
