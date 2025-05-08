# 6. 修饰根据修饰的分群获得Get top markers 修饰差异热图展示
# 函数：修饰根据修饰的分群获得Get top markers 修饰差异热图展示
# 参数：
#   m6A：单细胞修饰组数据
#   group：数据组名称
#   date：日期字符串
# 输出：
#   绘制并保存修饰差异热图，同时保存top markers数据
modification_markers_heatmap <- function(m6A, group, date) {
  library(RColorBrewer)
  label="m6Acluster"
  ncluster=length(unique(m6A$cluster_m6A))

  m6A$cluster_m6A <- factor(m6A$cluster_m6A, levels=paste0("sscm6A",1:ncluster))
  Idents(m6A) <- factor(m6A$cluster_m6A, levels=paste0("sscm6A",1:ncluster))
  m6A@active.ident <- factor(m6A$cluster_m6A, levels=paste0("sscm6A",1:ncluster))

  # top markers
  mk_m6A <- FindAllMarkers(m6A)
  markers <- mk_m6A
  top.markers <- markers %>%
               group_by(cluster) %>%
               dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.5) %>%
               top_n(5, wt=-log(p_val))

  # keep uniq gene
  top.markers <- top.markers[which(!top.markers$gene %in% names(table(top.markers$gene)[table(top.markers$gene) > 1])),] 

  # Reorder markers
  clusters_order_heatmap <- paste0("sscm6A",1:ncluster)
  top.markers$cluster <- factor(as.character(top.markers$cluster),levels=clusters_order_heatmap)
  top.markers <- top.markers[order(top.markers$cluster),]

  m6Acolors <- scales::hue_pal()(ncluster)
  names(m6Acolors) <- c(paste0("sscm6A",1:ncluster))

  CTcolors_heatmap         <- m6Acolors
  names(CTcolors_heatmap)  <- clusters_order_heatmap

  usedmarkers <- data.frame(rbind(top.markers, markers[which(markers$cluster=="sscm6A6"),]))
  usedmarkers <- top.markers %>%
               group_by(cluster) %>%
               top_n(5, wt=avg_log2FC)

  p <- DoHeatmap(object = m6A,features = rev(usedmarkers$gene),
                  assay = "peaks",slot="data",
                  group.colors=CTcolors,size = 3)+ scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#eb0520','#ef8a62','#fce1d2')), midpoint = 0, guide = "colourbar", aesthetics = "fill") +
                  theme(axis.text.y = element_text(face = "italic", size = 10),legend.position = "right")+
                  theme(# Hide panel borders and remove grid lines
                  #panel.border = element_blank(),
                  panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                  axis.line = element_line(colour = "black"))

  ggsave(filename = paste0("figs/",date,group,"-DoHeatmap_top5_",label,".pdf"),plot = p,width=15,height =15)
  saveRDS(usedmarkers, paste0("figs/",date,group,"-DoHeatmap_top5_",label,".rds"))
}

# 6. 修饰根据转录的分群获得Get top markers 差异热图展示
# 函数：修饰根据转录的分群获得Get top markers 差异热图展示
# 参数：
#   m6A_rna：包含转录组信息的单细胞修饰组数据
#   group：数据组名称
#   date：日期字符串
# 输出：
#   绘制并保存根据转录分群的修饰差异热图，同时保存top markers数据
transcription_markers_heatmap <- function(m6A_rna, group, date) {
  library(RColorBrewer)
  label="RNAcluster"
  ncluster=length(unique(na.omit(m6A_rna$cluster_RNA)))

  m6A_rna$cluster_RNA <- factor(m6A_rna$cluster_RNA, levels=paste0("sscRNA",1:ncluster))
  Idents(m6A_rna) <- factor(m6A_rna$cluster_RNA, levels=paste0("sscRNA",1:ncluster))
  m6A_rna@active.ident <- factor(m6A_rna$cluster_RNA, levels=paste0("sscRNA",1:ncluster))

  mk_rna <- FindAllMarkers(m6A_rna)
  markers <- mk_rna
  top.markers <- markers %>%
               group_by(cluster) %>%
               dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
               top_n(5, wt=-log(p_val))

  # keep uniq gene
  top.markers <- top.markers[which(!top.markers$gene %in% names(table(top.markers$gene)[table(top.markers$gene) > 1])),] 

  # Reorder markers
  clusters_order_heatmap <- paste0("sscRNA",1:ncluster)
  top.markers$cluster <- factor(as.character(top.markers$cluster),levels=clusters_order_heatmap)
  top.markers <- top.markers[order(top.markers$cluster),]

  # Reorder cells
  CTcolors <- RNAcolors
  names(CTcolors) <- paste0("sscRNA",1:ncluster)
  CTcolors_heatmap         <- CTcolors
  names(CTcolors_heatmap)  <- clusters_order_heatmap

  p <-
    DoHeatmap(m6A_rna, features = rev(top.markers$gene),
              assay = "peaks",
              slot = "data",
              group.colors=CTcolors,
              #angle = 45,
              group.by = c("cluster_RNA") ,
              label = T, disp.max = 2) &
    theme(axis.text.y = element_text(face = "italic", size = 10),legend.position = "right")

  ggsave(filename = paste0("figs/",date,group,"-DoHeatmap_top5_",label,".pdf"),plot = p,width=10,height =8)
  saveRDS(mk_rna, paste0("figs/",date,group,"-DoHeatmap_top5_",label,".rds"))
}
