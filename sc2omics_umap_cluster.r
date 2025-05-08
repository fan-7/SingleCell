###########################################
# 5. 单细胞转录组和单细胞修饰组合并
###########################################
# 设置工作目录，使用之前定义的工作目录变量workdir
setwd(file.path(workdir,"3_result"))
# 定义组合并的组名
group="IPInput"
# 创建组目录
dir.create(group)
# 切换到创建的组目录
setwd(group)
# 创建figs目录用于存储图表
dir.create("figs")

# 读取单细胞修饰组数据
m6A <- readRDS("IP_dedup-dimplot-umap-_TFIDF_LSI.rds")
# 读取单细胞转录组数据
rna <- readRDS("Input_dedup-dimplot-tsne-umap.rds")

# 向rna数据中添加m6A数据的元数据中的orig.ident和cluster_m6A字段
rna  <- AddMetaData(rna, m6A@meta.data[,c("orig.ident","cluster_m6A")])
# 向m6A数据中添加rna数据的元数据中的orig.ident和cluster_RNA字段
m6A  <- AddMetaData(m6A, rna@meta.data[,c("orig.ident","cluster_RNA")])

# 1. 修饰降维转录聚类
# 获取m6A数据中唯一的cluster_RNA的数量
ncluster=length(unique(na.omit(m6A$cluster_RNA)))
# 将m6A数据的标识设置为factor类型，以cluster_RNA为水平
Idents(m6A) <- factor(m6A$cluster_RNA, levels=paste0("sscRNA",1:ncluster))
# 将m6A数据的活动标识设置为factor类型，以cluster_RNA为水平
m6A@active.ident <- factor(m6A$cluster_RNA, levels=paste0("sscRNA",1:ncluster))

# 更新颜色，使用RColorBrewer包中的Set1调色板获取颜色，并命名
RNAcolors <- brewer.pal(ncluster, "Set1")
names(RNAcolors) <- c(paste0("sscRNA",1:ncluster))
# 修饰降维转录标签
label="m6AUMAP_RNAcluster"
# 绘制umap聚类图，设置点的大小、标题、颜色等
p_umap_cluster <- DimPlot(m6A, group.by="cluster_RNA",reduction = "umap", pt.size=0.5) +
  ggtitle("")+
  scale_color_manual(values = RNAcolors) +
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
# 保存umap聚类图为pdf文件
ggsave(filename = paste0("figs/",date,group,"-dimplot-umap_",label,".pdf"),plot = p_umap_cluster,width=6,height = 4)

# 2. 转录降维 转录聚类
# 从rna数据中提取与m6A数据相同细胞的子集
rna_m6A <- subset(rna, cells=colnames(m6A))
# 绘制umap聚类图，设置点的大小、标题、颜色等
p_umap_cluster <- DimPlot(rna_m6A, group.by="cluster_RNA", reduction = "umap", pt.size=0.5) +
  ggtitle("")+
  scale_color_manual(values = RNAcolors) +
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
# 保存umap聚类图为pdf文件
ggsave(filename = paste0("figs/",date,group,"-dimplot-umap_","_rnaUMAP_rnacluster",".pdf"),plot = p_umap_cluster,width=6,height = 4)

# 3. 转录降维 修饰聚类
# 获取rna_m6A数据中唯一的cluster_m6A的数量
ncluster=length(unique(rna_m6A$cluster_m6A))
# 将rna_m6A数据的标识设置为factor类型，以cluster_m6A为水平
Idents(rna_m6A) <- factor(rna_m6A$cluster_m6A, levels=paste0("sscm6A",1:ncluster))
# 将rna_m6A数据的活动标识设置为factor类型，以cluster_m6A为水平
rna_m6A@active.ident <- factor(rna_m6A$cluster_m6A, levels=paste0("sscm6A",1:ncluster))
# 将rna_m6A数据的orig.ident设置为RNA
rna_m6A$orig.ident <- "RNA"

# 获取颜色，使用scales包中的hue_pal函数获取颜色，并命名
m6Acolors <- scales::hue_pal()(ncluster)
names(m6Acolors) <- c(paste0("sscm6A",1:ncluster))

# 转录降维修饰标签
label="RNAUMAP_m6Acluster"
# 绘制umap聚类图，设置点的大小、标题、颜色等
p_umap_cluster <- DimPlot(rna_m6A, group.by="cluster_m6A", reduction = "umap", pt.size=0.5) +
  ggtitle("")+
  scale_color_manual(values = m6Acolors) +
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
        axis.line = element_line(colour = "black"))+
  theme(legend.direction = "horizontal")
# 保存umap聚类图为pdf文件
ggsave(filename = paste0("figs/",date,group,"-dimplot-umap_",label,".pdf"),plot = p_umap_cluster,width=10,height = 4)

# 4. 修饰降维 修饰聚类
# 获取颜色，使用scales包中的hue_pal函数获取颜色，并命名
m6Acolors <- scales::hue_pal()(ncluster)
names(m6Acolors) <- c(paste0("sscm6A",1:ncluster))
# 定义标签
label="_m6AUMAP_m6Acluster"
# 绘制umap聚类图，设置点的大小、标题、颜色等
p_umap_cluster <- DimPlot(m6A, reduction = "umap", group.by="cluster_m6A", pt.size=0.5) +
  ggtitle("")+
  scale_color_manual(values = m6Acolors) +
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
        axis.line = element_line(colour = "black"))+
  theme(legend.direction = "horizontal")
# 保存umap聚类图为pdf文件
ggsave(filename = paste0("figs/",date,group,"-dimplot-umap_",label,".pdf"),plot = p_umap_cluster,width=10,height = 4)
