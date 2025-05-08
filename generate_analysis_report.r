generate_analysis_report <- function(sobj_list) {
  # 创建结果目录
  result_dir <- file.path(workdir, "3_result")
  if(!dir.exists(result_dir)) dir.create(result_dir)
  
  # 保存关键结果
  saveRDS(sobj_list$clusters, 
          file.path(result_dir, "final_clusters.rds"))
  
}
