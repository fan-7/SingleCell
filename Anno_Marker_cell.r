# 1. marker注释热图
# 设置工作目录（假设workdir已定义）
workdir="IPbarcode/dedup"
setwd(workdir)
# 创建名为heatmap的目录
dir.create("heatmap")

# 定义组名
group='m6AUMAP_RNAcluster'
# 获取单细胞数据的count矩阵
matrix <- as.matrix(GetAssayData(sobj_wt_rna , slot = "counts"))
# 处理矩阵的行名
rownames(matrix) <- sapply(strsplit(rownames(matrix) , "\\."), "[[", 2)
# 运行自定义函数（假设run_sobj_fun已定义）
sobj <- run_sobj_fun(matrix)
# 设置分辨率
resolution=0.6
# 寻找所有标记
markers <- FindAllMarkers(sobj)
# 筛选显著的正标记
markers_pos <- markers[markers$p_val < 0.05 & markers$avg_log2FC > 0.5,]
# 将标记结果写入CSV文件
write.csv(x = markers, file = paste0("heatmap/", date,group,"-m6Acluster",resolution,"sobj_wt_gene_FindAllMarkers.csv"))

# 获取top标记
top.markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
  top_n(50, wt=-log(p_val))

# 保留唯一基因
top.markers <- top.markers[which(!top.markers$gene %in% names(table(top.markers$gene)[table(top.markers$gene) > 1])),] 

# 重新排序标记
clusters_order_heatmap <- c(paste0("sscRNA",1:11))
top.markers$cluster <- factor(as.character(top.markers$cluster),levels=clusters_order_heatmap)
top.markers <- top.markers[order(top.markers$cluster),]

# 重新排序细胞
Idents(sobj) <- factor(as.character(Idents(sobj)),levels=clusters_order_heatmap)

# 创建行注释数据框
row_annotation.df <- data.frame(markers=top.markers$cluster)
rownames(row_annotation.df) <- top.markers$gene

# 细胞聚类情况（假设sobj已定义）
col_annotation.df <- data.frame(clusters=Idents(sobj))
#col_annotation.df$GFP <- sobj$GFP
#col_annotation.df$Age <- sobj$Age

# 颜色定义
CTcolors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99")
names(CTcolors) <- c(paste0("sscRNA",1:11))

# 绘制热图
p <-
  DoHeatmap(m6A, features = rev(top.markers$gene),
            assay = "RNA",
            group.colors=CTcolors,
            #angle = 45,
            group.by = c("seurat_clusters") ,
            #group.colors = Fibs_col.palette,
            label = T, disp.max = 2) &
  theme(axis.text.y = element_text(face = "italic", size = 10),legend.position = "right")

# 保存热图为PDF文件
ggsave(filename = paste0("heatmap/",date,group,"-DoHeatmap_top100-keepdup_",filter,label,".pdf"),plot = p,width=4,height = 4)

# 1. 对差异的marker进行注释和特征表达 (选择差异修饰的marker)
# 筛选含有指定修饰的基因（假设gene_to_keep已定义）
#gene_to_keep <- "Exoc4"
#cells_with_gene <- rownames(sobj_wt_rna@assays$RNA@counts)[grep(gene_to_keep, rownames(sobj_wt_rna@assays$RNA@counts))]

# 在 FeaturePlot 中可视化筛选后的细胞
grep("Gm22068", rownames(sobj_wt_rna))
p <- FeaturePlot(sobj_wt_rna,"ENSMUSG00000096205.Gm22068.snRNA", reduction="tsne", cols = c('#CCCCCC','#FF0000'),pt.size=0.5,min.cutoff = 1) + ggtitle("ENSMUSG00000096205.Gm22068.snRNA") + NoAxes()
# 保存特征图为PDF文件
ggsave(plot = p,filename = paste0("wt/ENSMUSG00000096205.Gm22068.snRNA.pdf"),width=4,height=4)

# 20230912 展示单细胞多组学细胞对应的结果
# 切换工作目录（假设路径已去除个人信息）
setwd("hitwtkopeak/sc/dedup/wt")
# 创建名为matchcell的目录
dir.create("matchcell")
# 读取单细胞修饰组数据
m6A <- readRDS("sce231008/E10_5-IP/qc_tag/star/count/hitwtkopeak/sc/dedup/wt/October-13-23E10_5-ZLZ-dimplot-tsne-dedup_peakres06umap.rds")
# 读取单细胞转录组数据
rna <- readRDS("sce231008/E10_5-ZLZ/qc_tag/star/count/IPbarcode/dedup/figs/October-13-23E10_5-ZLZ-dimplot-tsne-keepdup_dedupgeneres0.6.rds")
# 给转录组数据添加模态信息
rna$modality <- "RNA"
# 给修饰组数据添加模态信息
m6A$modality <- "m6A"
# 将转录组和修饰组数据组合成列表
seurat <- list(rna, m6A)

# 定义模态数量
nmodalities      = 2
# 初始化umap嵌入列表
umap_embeddings <- list()

# 提取并处理每个模态的umap嵌入
for(i in 1:length(seurat[1:nmodalities])){
  umap_embeddings[[i]]               <- as.data.frame(seurat[[i]]@reductions[['tsne']]@cell.embeddings)
  umap_embeddings[[i]]$tSNE_1        <- umap_embeddings[[i]]$tSNE_1 + (i-1)*60
  umap_embeddings[[i]]$modality      <- unique(seurat[[i]]$modality)
  umap_embeddings[[i]]$cluster       <- seurat[[i]]$RNA_snn_res.0.6
  umap_embeddings[[i]]$cell_barcode  <- rownames(umap_embeddings[[i]])
}

# 合并umap嵌入数据
umap.embeddings.merge               <- purrr::reduce(umap_embeddings,rbind)
# 找到共同的细胞
common.cells                        <- table(umap.embeddings.merge$cell_barcode)
common.cells                        <- names(common.cells[common.cells==nmodalities])

# 筛选出共同细胞的数据
umap.embeddings.merge               <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]

# 计算细胞数量
ncell=length(common.cells)
# 随机选择细胞
cells.connect <- sample(umap.embeddings.merge$cell_barcode,ncell)

# 绘制带有连接的umap图
p1 <- ggplot(data=umap.embeddings.merge,aes(x=tSNE_1,y=tSNE_2,col=cluster)) +
  geom_point(size=0.2) +
  scale_color_brewer(palette = 'Set1') +
  geom_line(data=umap.embeddings.merge, aes(group=cell_barcode,col=cluster),alpha=0.1,size=0.01) +
  theme_classic() + NoAxes() + NoLegend()

# 保存umap图为PDF文件
ggsave(filename = paste0('matchcell/tsne_with_connections_5x5.pdf'), plot=p1, width=5,height=5,bg='transparent')

# 20231013 转录组注释细胞类型
# 切换工作目录（假设路径已去除个人信息）
setwd("202310result/sce231008E105/Input")
# 读取CSV文件
input <- read.csv("October-13-23peak_gene-m6Acluster0.2sobj_FindAllMarkers.csv")
# 处理输入数据的mk列
input$mk <- sapply(strsplit(input$gene , "\\."), "[[", 2)


################################################################
# 定义工作目录（去除个人路径，假设后续会正确设置）
# workdir="sce231008/E10_5-ZLZ/qc_tag/star/count/IPbarcode/dedup"
# setwd(workdir)
# 创建 heatmap 目录
dir.create("heatmap")

# 定义组名
group='m6AUMAP_RNAcluster'

# 获取数据矩阵并处理行名
matrix <- as.matrix(GetAssayData(sobj_wt_rna , slot = "counts"))
rownames(matrix) <- sapply(strsplit(rownames(matrix) , "\\."), "[[", 2)

# 运行自定义函数（假设 run_sobj_fun 已定义）
sobj <- run_sobj_fun(matrix)

# 设置分辨率
resolution=0.6

# 寻找所有标记
markers <- FindAllMarkers(sobj)

# 筛选显著的正标记
markers_pos <- markers[markers$p_val < 0.05 & markers$avg_log2FC > 0.5,]

# 将标记结果写入CSV文件
write.csv(x = markers, file = paste0("heatmap/", date,group,"-m6Acluster",resolution,"sobj_wt_gene_FindAllMarkers.csv"))

# 获取 top 标记
top.markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
  top_n(50, wt=-log(p_val))

# 保留唯一基因
top.markers <- top.markers[which(!top.markers$gene %in% names(table(top.markers$gene)[table(top.markers$gene) > 1])),] 

# 重新排序标记
clusters_order_heatmap <- c(paste0("sscRNA",1:11))
top.markers$cluster <- factor(as.character(top.markers$cluster),levels=clusters_order_heatmap)
top.markers <- top.markers[order(top.markers$cluster),]

# 重新排序细胞
Idents(sobj) <- factor(as.character(Idents(sobj)),levels=clusters_order_heatmap)

# 创建行注释数据框
row_annotation.df <- data.frame(markers=top.markers$cluster)
rownames(row_annotation.df) <- top.markers$gene

# 细胞聚类情况
col_annotation.df <- data.frame(clusters=Idents(sobj))
#col_annotation.df$GFP <- sobj$GFP
#col_annotation.df$Age <- sobj$Age

# 颜色定义
CTcolors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F","#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99")
names(CTcolors) <- c(paste0("sscRNA",1:11))

# 绘制热图
p <-
  DoHeatmap(m6A, features = rev(top.markers$gene),
            assay = "RNA",
            group.colors=CTcolors,
            #angle = 45,
            group.by = c("seurat_clusters") ,
            #group.colors = Fibs_col.palette,
            label = T, disp.max = 2) &
  theme(axis.text.y = element_text(face = "italic", size = 10),legend.position = "right")

# 保存热图为PDF文件
ggsave(filename = paste0("heatmap/",date,group,"-DoHeatmap_top100-keepdup_",filter,label,".pdf"),plot = p,width=4,height = 4)

# 在 FeaturePlot 中可视化筛选后的细胞
grep("Gm22068", rownames(sobj_wt_rna))
p <- FeaturePlot(sobj_wt_rna,"ENSMUSG00000096205.Gm22068.snRNA", reduction="tsne", cols = c('#CCCCCC','#FF0000'),pt.size=0.5,min.cutoff = 1) + ggtitle("ENSMUSG00000096205.Gm22068.snRNA") + NoAxes()
ggsave(plot = p,filename = paste0("wt/ENSMUSG00000096205.Gm22068.snRNA.pdf"),width=4,height=4)

# 20230912 展示单细胞多组学细胞对应的结果
# 设置工作目录（去除个人路径，假设后续会正确设置）
# setwd("sce231008/E10_5-IP/qc_tag/star/count/hitwtkopeak/sc/dedup/wt")
# 创建 matchcell 目录
dir.create("matchcell")

# 读取单细胞修饰组数据
m6A <- readRDS("October-13-23E10_5-ZLZ-dimplot-tsne-dedup_peakres06umap.rds")
# 读取单细胞转录组数据
rna <- readRDS("October-13-23E10_5-ZLZ-dimplot-tsne-keepdup_dedupgeneres0.6.rds")

# 添加模态信息
rna$modality <- "RNA"
m6A$modality <- "m6A"

# 将转录组和修饰组数据组合成列表
seurat <- list(rna, m6A)

# 定义模态数量
nmodalities = 2
# 初始化umap嵌入列表
umap_embeddings <- list()

# 提取并处理每个模态的umap嵌入
for(i in 1:length(seurat[1:nmodalities])){
  umap_embeddings[[i]] <- as.data.frame(seurat[[i]]@reductions[['tsne']]@cell.embeddings)
  umap_embeddings[[i]]$tSNE_1 <- umap_embeddings[[i]]$tSNE_1 + (i-1)*60
  umap_embeddings[[i]]$modality <- unique(seurat[[i]]$modality)
  umap_embeddings[[i]]$cluster <- seurat[[i]]$RNA_snn_res.0.6
  umap_embeddings[[i]]$cell_barcode <- rownames(umap_embeddings[[i]])
}

# 合并umap嵌入数据
umap.embeddings.merge <- purrr::reduce(umap_embeddings,rbind)

# 找到共同的细胞
common.cells <- table(umap.embeddings.merge$cell_barcode)
common.cells <- names(common.cells[common.cells==nmodalities])

# 筛选出共同细胞的数据
umap.embeddings.merge <- umap.embeddings.merge[umap.embeddings.merge$cell_barcode %in% common.cells,]

# 计算细胞数量
ncell=length(common.cells)
# 随机选择细胞
cells.connect <- sample(umap.embeddings.merge$cell_barcode,ncell)

# 绘制带有连接的umap图
p1 <- ggplot(data=umap.embeddings.merge,aes(x=tSNE_1,y=tSNE_2,col=cluster)) +
  geom_point(size=0.2) +
  scale_color_brewer(palette = 'Set1') +
  geom_line(data=umap.embeddings.merge, aes(group=cell_barcode,col=cluster),alpha=0.1,size=0.01) +
  theme_classic() + NoAxes() + NoLegend()

# 保存umap图为PDF文件
ggsave(filename = paste0('matchcell/tsne_with_connections_5x5.pdf'), plot=p1, width=5,height=5,bg='transparent')

# 20231013 转录组注释细胞类型
# 设置工作目录（去除个人路径，假设后续会正确设置）
# setwd("Input")
# 读取CSV文件
input <- read.csv("October-13-23peak_gene-m6Acluster0.2sobj_FindAllMarkers.csv")

# 处理输入数据的mk列
input$mk <- sapply(strsplit(input$gene , "\\."), "[[", 2)

# 读取细胞类型标记文件
nature <- read.table("2019nature_embryo_marker.txt", sep="\t", header = T)
# 提取自然标记
nature_mk <- unlist(strsplit(nature$Markers.used.for.cell.type.identification, ", "))
# 筛选input中与自然标记匹配的行
input[which(input$mk %in% nature_mk), ]

# 读取细胞类型标记文件
nature <- read.table("../../../202310paper/2019nature_embryo_marker.txt", sep="\t", header = T)
# 提取自然标记
nature_mk <- unlist(strsplit(nature$Markers.used.for.cell.type.identification, ", "))
# 筛选input中与自然标记匹配的行
input[which(input$mk %in% nature_mk), ]
