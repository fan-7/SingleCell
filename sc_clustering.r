sc_clustering <- function(sobj, res=0.6) {
  # 标准化流程
  sobj <- NormalizeData(sobj) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
  
  # 非线性降维
  sobj <- RunUMAP(sobj, dims=1:20) %>%
    RunTSNE()
  
  # 聚类分析
  sobj <- FindNeighbors(sobj, dims=1:20) %>%
    FindClusters(resolution=res)
  
  return(sobj)
}
