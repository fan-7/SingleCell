###########################################################
# 1. 数据复制模块（假设在 Shell 环境下）
###########################################################
#!/bin/bash
# 目的：将数据从源位置复制到目标位置
# 定义目标路径
target_dir="/shared_data/sscm6A/0_data/"

# 复制数据
# 示例数据复制，可根据实际情况修改源路径和文件名
cp /source_path/keepdup_sobj1 $target_dir/input_keepdup.rds
cp /source_path/dedup_sobj1 $target_dir/input_dedup.rds
cp /source_path/IP_dedup_sobj1 $target_dir/IP_dedup.rds
cp /source_path/IP_keepdup_sobj1 $target_dir/IP_keepdup.rds
###########################################################
#2. R 环境设置与包加载模块
###########################################################
# 目的：设置R环境并加载所需的R包
# 说明：创建工作目录并加载单细胞分析常用包

# 设置工作目录，实际使用时修改为合适路径
workdir="/shared_data/sscm6A/"
setwd(workdir)

# 加载所需的R包
library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(tidyverse)
library(CytoTRACE)
library(Signac)
library(viridis)
library(org.Mm.eg.db)
library(pheatmap)
library(ggridges)
library(gridExtra)
###########################################################
#3. 单细胞捕获情况分析模块
###########################################################
# 目的：分析单细胞捕获情况，包括每个细胞的UMI和基因数量分布
# 输入：group - 数据组名称，如"Input_keepdup"
# 输出：生成相关可视化图表并保存

analyze_capture <- function(group) {
  date=format(Sys.time(), format="%B-%d-%y")
  sobj <- readRDS(paste0(workdir, "/0_data/",group,".rds"))
  setwd(file.path(workdir,"3_result"))
  dir.create(group)
  setwd(group)
  
  sobj$UMIs = sobj$nCount_RNA
  sobj$Genes = sobj$nFeature_RNA
  me2 <- sobj@meta.data %>%
    mutate(UmisPerGene = UMIs/Genes) %>%
    pivot_longer(cols = c(UMIs, Genes), names_to = "Type", values_to = "Number") %>% data.frame()
  
  HMcolors <- c("#D6604D", "#4393C3")
  names(HMcolors) <- c("UMIs","Genes")
  dir.create("UMIvsGene")
  setwd("UMIvsGene/")
  dir.create("figures")
  
  # 每个细胞的密度图
  plot_density <- me2 %>% ggplot(aes(x=Number,fill =Type)) +
    geom_density(alpha = 0.7) +
    scale_x_log10() +
    scale_color_manual(values = HMcolors)  +
    scale_fill_manual(values = HMcolors) +
    ggthemes::theme_few() +
    ylab("Cell density") +
    theme(axis.text.x = element_text(size=15,color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title.x =element_text(size=15, color = "black"),
          axis.title.y =element_text(size=15, color = "black"),
          axis.ticks.length.x=unit(0.2, "cm"),
          axis.ticks.length.y=unit(0.2, "cm"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 8),
          axis.ticks.y=element_line(color="black",size=1,lineend = 8)) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  ggsave(paste0("./figures/densityplot_ReadsvsGene_perCell.",group,".pdf"),plot = plot_density,width=4,height = 4)
  
  # 每个细胞基因箱线图
  plot_box <- me2 %>%
    ggplot(aes(x=Type, y=log10(Number), fill=Type)) +
    geom_violin(aes(fill=Type)) +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    ggthemes::theme_few() +
    scale_color_manual(values = HMcolors)  +
    scale_fill_manual(values = HMcolors)  +
    theme(axis.text.x = element_text(size=15,color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title.x =element_text(size=15, color = "black"),
          axis.title.y =element_text(size=15, color = "black"),
          axis.ticks.length.x=unit(0.2, "cm"),
          axis.ticks.length.y=unit(0.2, "cm"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 8),
          axis.ticks.y=element_line(color="black",size=1,lineend = 8)) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  ggsave(paste0("./figures/boxplot_ReadsvsGene_perCell.",group,".pdf"),plot = plot_box,width=4,height = 4)
  
  saveRDS(me2, paste0("./figures/boxplot_ReadsvsGene_perCell.",group,".rds"))
}
###########################################################
#4. 单细胞聚类分析模块
###########################################################
# 目的：进行单细胞聚类分析
# 输入：group - 数据组名称，如"Input_dedup"
# 输出：生成聚类可视化图表并保存

cluster_single_cell <- function(group) {
  date=format(Sys.time(), format="%B-%d-%y")
  sobj <- readRDS(paste0(workdir, "/0_data/",group,".rds"))
  setwd(file.path(workdir,"3_result"))
  dir.create(group)
  setwd(group)
  
  # 定义聚类函数
  run_sobj_fun <- function(count_matrix){
    RNA.sobj <- CreateSeuratObject(counts = count_matrix, assay = 'RNA')
    all.genes <- rownames(RNA.sobj)
    #1 标准化数据， 对数标准化 尺度因子是1万
    RNA.sobj <- NormalizeData(RNA.sobj, normalization.method = "LogNormalize", scale.factor = 10000)
    #2 选择特征变量，通过vst方法选择2000特征变量
    RNA.sobj <- FindVariableFeatures(RNA.sobj, selection.method = "vst", nfeatures = 2000)
    #3 归一化数据
    RNA.sobj <- ScaleData(RNA.sobj, features = all.genes)
    #4 PCA降维 基于2000的特征计算
    RNA.sobj <- RunPCA(RNA.sobj, features = VariableFeatures(object = RNA.sobj))
    #5 UMAP聚类 维度20
    RNA.sobj <- RunUMAP(RNA.sobj, dims = 1:20)
    #6 TSNE聚类
    RNA.sobj <-  RunTSNE(RNA.sobj, reduction = "pca", seed.use = 1,tsne.method = "Rtsne", check_duplicates = FALSE)
    #7  K最近邻
    RNA.sobj <- FindNeighbors(object = RNA.sobj,
                              dims = 1:20)
    #8. 确定聚类 resolution
    RNA.sobj <- FindClusters(object = RNA.sobj,
                             resolution = 0.6)
    return(RNA.sobj)
  }
  
  cm <- GetAssayData(sobj, slot="counts")
  sobj_ct<- run_sobj_fun(cm)
  dir.create("figs")
  sobj_ct$orig.ident <- "rna"
  sobj_ct$cluster_RNA <- paste0("sscRNA",as.numeric(as.character(sobj_ct$seurat_clusters))+1)
  
  # TSNE聚类可视化
  p_tsne_cluster <- DimPlot(sobj_ct, reduction = "tsne", group.by = "cluster_RNA", pt.size=0.5) +
    ggtitle("")+
    scale_color_brewer(palette = 'Set1') +
    ggthemes::theme_few() +
    theme(axis.text.x = element_text(size=15,color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title.x =element_text(size=15, color = "black"),
          axis.title.y =element_text(size=15, color = "black"),
          axis.ticks.length.x=unit(0.2, "cm"),
          axis.ticks.length.y=unit(0.2, "cm"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 8),
          axis.ticks.y=element_line(color="black",size=1,lineend = 8)) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  ggsave(filename = paste0("figs/",date,group,"-dimplot-tsne-",".pdf"),plot = p_tsne_cluster,width=6,height = 4)
  
  # UMAP聚类可视化
  p_umap_cluster <- DimPlot(sobj_ct, reduction = "umap", group.by = "cluster_RNA", pt.size=0.5) +
    ggtitle("")+
    scale_color_brewer(palette = 'Set1') +
    ggthemes::theme_few() +
    theme(axis.text.x = element_text(size=15,color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title.x =element_text(size=15, color = "black"),
          axis.title.y =element_text(size=15, color = "black"),
          axis.ticks.length.x=unit(0.2, "cm"),
          axis.ticks.length.y=unit(0.2, "cm"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 8),
          axis.ticks.y=element_line(color="black",size=1,lineend = 8)) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  ggsave(filename = paste0("figs/",date,group,"-dimplot-umap-",".pdf"),plot = p_umap_cluster,width=6,height = 4)
  saveRDS(sobj_ct, paste0("figs/",date,group,"-dimplot-tsne-umap",".rds"))
  
  # 使用scSHC确定显著分类的结果
  shc <- scSHC(GetAssayData(sobj_ct, slot="data"), batch = NULL, alpha = 0.05, num_features = 2000, num_PCs = 50, parallel = T, cores = 12)
  sobj_shc <- AddMetaData(sobj_ct, data.frame(shc[[1]]))
  
  p_umap_cluster_shc <- DimPlot(sobj_shc, group.by="shc..1..", reduction = "umap",pt.size=0.1) +
    ggtitle("")+
    scale_color_brewer(palette = 'Set1')+
    ggthemes::theme_few() +
    theme(axis.text.x = element_text(size=15,color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title.x =element_text(size=15, color = "black"),
          axis.title.y =element_text(size=15, color = "black"),
          axis.ticks.length.x=unit(0.2, "cm"),
          axis.ticks.length.y=unit(0.2, "cm"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 8),
          axis.ticks.y=element_line(color="black",size=1,lineend = 8)) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  ggsave(filename = paste0("figs/",date,group,"-dimplot-umap-","_shc005",".pdf"),plot = p_umap_cluster_shc,width=4,height = 4)
  saveRDS(sobj_shc, paste0("figs/",date,group,"-dimplot-umap-","_shc005",".rds"))
}
###########################################################
#5. 单细胞修饰组降维展示模块
###########################################################
# 目的：对单细胞修饰组进行降维展示
# 输入：group - 数据组名称，如"IP_dedup"
# 输出：生成修饰组降维可视化图表并保存

visualize_modification <- function(group) {
  date=format(Sys.time(), format="%B-%d-%y")
  sobj <- readRDS(paste0(workdir, "/0_data/",group,".rds"))
  setwd(file.path(workdir,"3_result"))
  dir.create(group)
  setwd(group)
  dir.create("figs")
  
  counts <- GetAssayData(sobj, slot="counts")
  metadata <- sobj@meta.data
  sobj_peak <- CreateSeuratObject(
    counts = counts,
    assay = 'peaks',
    project = 'sscm6A',
    min.cells = 1,
    meta.data = metadata
  )
  
  sobj_peak <- RunTFIDF(sobj_peak)
  sobj_peak <- FindTopFeatures(sobj_peak, min.cutoff = 'q0')
  sobj_peak <- RunSVD(
    object = sobj_peak,
    assay = 'peaks',
    reduction.key = 'LSI_',
    reduction.name = 'lsi'
  )
  
  sobj_peak <- RunUMAP(object = sobj_peak, reduction = 'lsi', dims = 1:30)
  sobj_peak <- FindNeighbors(object = sobj_peak, reduction = 'lsi', dims = 1:30)
  sobj_peak <- FindClusters(object = sobj_peak, verbose = FALSE, algorithm = 3)
  
  sobj_peak$orig.ident <- "m6A"
  sobj_peak$cluster_m6A <- paste0("sscm6A",as.numeric(as.character(sobj_peak$seurat_clusters))+1)
  
  CTcolors <- scales::hue_pal()(length(unique(sobj_peak$seurat_clusters)))
  names(CTcolors) <- paste0("sscm6A",1:(length(unique(sobj_peak$seurat_clusters))))
  
  p_umap_cluster <- DimPlot(sobj_peak, group.by="cluster_m6A", reduction = "umap",pt.size=0.5) +
    ggtitle("")+
    scale_color_manual(values = CTcolors) +
    ggthemes::theme_few() +
    theme(axis.text.x = element_text(size=15,color = "black"),
          axis.text.y = element_text(size=15, color = "black"),
          axis.title.x =element_text(size=15, color = "black"),
          axis.title.y =element_text(size=15, color = "black"),
          axis.ticks.length.x=unit(0.2, "cm"),
          axis.ticks.length.y=unit(0.2, "cm"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 8),
          axis.ticks.y=element_line(color="black",size=1,lineend = 8)) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype
