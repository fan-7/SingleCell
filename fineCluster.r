# 7. 单细胞修饰精细分群
# 这个脚本主要进行单细胞修饰精细分群，包括转录类群的细分和修饰类群细分的转录，并绘制相应的 umap 图进行可视化，最后在 FeaturePlot 中可视化筛选后的细胞。在使用时，请确保相关的数据对象（如rna_m6A, m6A_rna, sobj_wt）和变量（如m6Acolors, date, group）已正确定义。
# 设置工作目录（假设work_dir已定义）
setwd(work_dir)
# 创建名为fineCluster的目录
dir.create("fineCluster")

# 转录类群的细分
# 设置颜色（假设m6Acolors已定义）
CTcolors <- m6Acolors
label="_RNAUMAP_m6Acluster"
# 假设rna_m6A为单细胞数据对象
# 将rna2列赋值为"Others"
rna_m6A$rna2 <- "Others"
# 定义新的类群向量
newCT <-  c("sscm6A8", "sscm6A12", "sscm6A13", "sscm6A16", "sscm6A17","sscm6A19")
# 找到满足条件的索引
index <- which((rna_m6A$cluster_m6A %in% newCT) & (rna_m6A$cluster_RNA %in% c("sscRNA2")))
# 根据索引更新rna2列的值
rna_m6A$rna2[index] <- rna_m6A$cluster_m6A[index]
# 从m6Acolors中提取对应新类群的颜色
subcolors <- m6Acolors[newCT]
# 命名颜色
names(subcolors) <- newCT

# 绘制umap图
p_umap_cluster <- DimPlot(rna_m6A, reduction = "umap", group.by="rna2", pt.size=0.2) +
  ggtitle("")+
  scale_color_manual(values = subcolors) +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15, color = "black"),
        axis.title.x =element_text( size=15, color = "black"),
        axis.title.y =element_text(size=15, color = "black"),
        axis.ticks.length.x=unit(0.2, "cm"),
        axis.ticks.length.y=unit(0.2, "cm"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 8),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 8)) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

# 保存umap图为pdf文件
ggsave(filename = paste0("fineCluster/",date,group,"-dimplot-umap_RNA2_",label,".pdf"),plot = p_umap_cluster,width=6,height = 4)
# 保存数据对象
saveRDS(rna_m6A, paste0("fineCluster/",date,group,"-dimplot-umap_",label,".rds"))

# 修饰类群细分的转录
label="_m6AUMAP_RNAcluster"
# 假设m6A_rna为单细胞数据对象
# 将rna2列赋值为"Others"
m6A_rna$rna2 <- "Others"
# 找到满足条件的索引
index <- which((m6A_rna$cluster_m6A %in% newCT) & (m6A_rna$cluster_RNA %in% c("sscRNA2")))
# 根据索引更新rna2列的值
m6A_rna$rna2[index] <- as.character(m6A_rna$cluster_RNA)[index]
# 定义新的颜色向量
subcolors <- c("#377EB8", "gray")
# 命名颜色
names(subcolors) <- c("sscRNA2", "Others")

# 绘制umap图
p_umap_cluster <- DimPlot(m6A_rna, reduction = "umap", group.by="rna2", pt.size=0.5) +
  ggtitle("")+
  scale_color_manual(values = subcolors) +
  ggthemes::theme_few() +
  theme(axis.text.x = element_text(size=15,color = "black"),
        axis.text.y = element_text(size=15, color = "black"),
        axis.title.x =element_text( size=15, color = "black"),
        axis.title.y =element_text(size=15, color = "black"),
        axis.ticks.length.x=unit(0.2, "cm"),
        axis.ticks.length.y=unit(0.2, "cm"),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 8),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 8)) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

# 保存umap图为pdf文件
ggsave(filename = paste0("fineCluster/",date,group,"-dimplot-umap_",label,".pdf"),plot = p_umap_cluster,width=6,height = 4)
# 保存数据对象
saveRDS(m6A_rna, paste0("fineCluster/",date,group,"-dimplot-umap_",label,".rds"))

# 在 FeaturePlot 中可视化筛选后的细胞
# 假设sobj_wt为单细胞数据对象
# 找到含有"Gm22068"的行名
grep("Gm22068", rownames(sobj_wt_rna))
# 绘制特征图
p <- FeaturePlot(sobj_wt,"peak.11.87426739.87426945.Gm22068.snRNA", reduction="tsne", cols = c('#CCCCCC','#FF0000'),pt.size=0.5,min.cutoff = 1) + ggtitle("peak.11.87426739.87426945.Gm22068.snRNA") + NoAxes()
# 保存特征图为pdf文件
ggsave(plot = p,filename = paste0("wt/peak.11.87426739.87426945.Gm22068.snRNA.pdf"),width=4,height=4)
