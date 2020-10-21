library(Seurat)
BiocManager::install('Seurat')
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
setwd("D://data/ScRNAcode/wang/")
##=======================1.创建Seurat对象========================
dir <- 'filtered_gene_bc_matrices/ref/'
counts <- Read10X(dir)
wang = CreateSeuratObject(counts, project = "zxz", min.cells=3, min.features = 200)
dim(wang)

##=======================2.数据质控与标准化================================
##dir.create('QC')
##提取线粒体基因
wang[["percent.mt"]] <- PercentageFeatureSet(wang, pattern='^ATMG')
violin <- VlnPlot(wang,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 3)
ggsave("QC/vlnplot-before-qc.pdf", plot = violin, width = 15, height = 6) 
ggsave("QC/vlnplot-before-qc.png", plot = violin, width = 15, height = 6) 
plot1 <- FeatureScatter(wang, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("QC/pearplot-before-qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot-before-qc.png", plot = pearplot, width = 12, height = 5)
##设置质控标准
wang<-subset(wang,subset=nFeature_RNA>500 & nFeature_RNA<5000 &percent.mt<0.5)
dim(wang)
## 绘制质量控制后的图
violin <-VlnPlot(wang,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 pt.size = 0.1, 
                 ncol = 3)
ggsave("QC/vlnplot-after-qc.pdf", plot = violin, width = 15, height = 6) 
ggsave("QC/vlnplot-after-qc.png", plot = violin, width = 15, height = 6)
## 基因表达量标准化
## 它的作用是让测序数据量不同的细胞的基因表达量具有可比性。计算公式如下：
## 标准化后基因表达量 = log1p（10000*基因counts/细胞总counts）
wang <- NormalizeData(wang, normalization.method = "LogNormalize", scale.factor = 10000)
##=======================3.数据降维与聚类==================================
## 寻找高变基因
## dir.create("cluster")
wang <- FindVariableFeatures(wang,mean.cutoff=c(0.0125,3),dispersion.cutoff =c(1.5,Inf) )
top10 <- head(VariableFeatures(wang), 10)
plot1 <- VariableFeaturePlot(wang) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
## 横坐标是某基因在所有细胞中的平均表达值，纵坐标是此基因的方差。
## 红色的点是被选中的高变基因，黑色的点是未被选中的基因，变异程度最高的10个基因在如图中标注了基因名称。
ggsave("cluster/VariableFeatures.pdf", plot = plot, width = 8, height = 6) 
ggsave("cluster/VariableFeatures.png", plot = plot, width = 8, height = 6)
## 数据缩放
scale.genes <-  rownames(wang)
wang <- ScaleData(wang, features = scale.genes)
## PCA降维并提取主成分
wang <- RunPCA(wang, features = VariableFeatures(wang),npcs = 100) 
plot1 <- DimPlot(wang, reduction = "pca") 
plot2 <- ElbowPlot(wang, ndims=40, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster/pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("cluster/pca.png", plot = plotc, width = 8, height = 4)
## 细胞聚类
## 此步利用 细胞-PC值 矩阵计算细胞之间的距离，
## 然后利用距离矩阵来聚类。其中有两个参数需要人工选择，
## 第一个是FindNeighbors()函数中的dims参数，需要指定哪些pc轴用于分析，选择依据是之前介绍的cluster/pca.png文件中的右图。
## 第二个是FindClusters()函数中的resolution参数，需要指定0.1-1.0之间的一个数值，用于决定clusters的相对数量，数值越大cluters越多。
wang <- FindNeighbors(object = wang, dims = 1:100)
wang <- FindClusters(object = wang, resolution = 1.0)
table(wang@meta.data$seurat_clusters)
## 非线性降维
## tsne
wang <- RunTSNE(wang, dims =1:40)
embed_tsne <- Embeddings(wang, 'tsne')
write.csv(embed_tsne,'cluster/embed_tsne.csv')
plot1 = DimPlot(wang, reduction = "tsne" ,label = "T", pt.size = 1,label.size = 4)
ggsave("cluster/tSNE_cluster.pdf", plot = plot1, width = 8, height = 7)
ggsave("cluster/tSNE_cluster.png", plot = plot1, width = 8, height = 7)

## UMAP
wang <- RunUMAP(wang,n.neighbors = 30,metric = 'correlation',min.dist = 0.3,dims = 1:40)
embed_umap <- Embeddings(wang, 'umap')
write.csv(embed_umap,'cluster/embed_umap.csv')
plot2 = DimPlot(wang, reduction = "umap",label = "T", pt.size = 1,label.size = 4) 
ggsave("cluster/UMAP_cluster.pdf", plot = plot2, width = 8, height = 7)
ggsave("cluster/UMAP_cluster.png", plot = plot2, width = 8, height = 7)
## 合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("cluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("cluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)


##========================4.特异基因作图==============================
## dir.create("cell_identify")
pdf("cell_identify/specific_vlnplot.pdf",height = 14,width = 14)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G78370')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G71140')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G37160')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G56870')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G45700')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G22080')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G14020')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G11550')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G28790')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G73590')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G28500')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G26760')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G44030')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G23410')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G27140')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G21750')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G04025')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G04890')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G14650')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G17280')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G29025')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G61470')
CombinePlots(plots = list(plot1,plot2),legend = 'none',ncol=1)
dev.off()
pdf("cell_identify/fig1_vlnplot.pdf",height = 14,width = 14)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G22710')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G25710')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G79580')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G49270')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G57620')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G63090')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G13620')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G20840')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G37090')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G42840')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G77690')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G54890')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
dev.off()
VlnPlot(object = wang,pt.size=0.1, features ='AT1G14960')

##==============================5.修改聚类标号=====================
##修改聚类号重新做图
new.cluster.ids<-c("2",'1','4','5','13','3','12','21','8','6','11',
                   '9','7','10','6','15','22','14','17','19','16',
                   '20','18','23','24')
names(new.cluster.ids) <- levels(wang)
wang <- RenameIdents(wang, new.cluster.ids)
Idents(wang)<-factor(Idents(wang),levels=mixedsort(levels(Idents(wang))))
wang <- RunTSNE(wang, dims =1:40)
embed_tsne <- Embeddings(wang, 'tsne')
write.csv(embed_tsne,'cluster/embed_tsne-new.csv')
plot1 = DimPlot(wang, reduction = "tsne" ,label = "T", pt.size = 1,label.size = 4)
ggsave("cluster/tSNE_cluster-new.pdf", plot = plot1, width = 8, height = 7)
ggsave("cluster/tSNE_cluster-new.png", plot = plot1, width = 8, height = 7)

## UMAP
wang <- RunUMAP(wang,n.neighbors = 30,metric = 'correlation',min.dist = 0.3,dims = 1:40)
embed_umap <- Embeddings(wang, 'umap')
write.csv(embed_umap,'cluster/embed_umap-new.csv')
plot2 = DimPlot(wang, reduction = "umap",label = "T", pt.size = 1,label.size = 4) 
ggsave("cluster/UMAP_cluster.pdf", plot = plot2, width = 8, height = 7)
ggsave("cluster/UMAP_cluster.png", plot = plot2, width = 8, height = 7)
## 合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("cluster/tSNE_UMAP-new.pdf", plot = plotc, width = 10, height = 5)
ggsave("cluster/tSNE_UMAP-new.png", plot = plotc, width = 10, height = 5)

pdf("cell_identify/specific_vlnplot-new.pdf",height = 14,width = 14)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G78370')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G71140')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G37160')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G56870')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G45700')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G22080')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G14020')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G11550')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G28790')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G73590')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G28500')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G26760')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G44030')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G23410')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G27140')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G21750')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G04025')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G04890')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G14650')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT4G17280')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G29025')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G61470')
CombinePlots(plots = list(plot1,plot2),legend = 'none',ncol=1)
dev.off()
pdf("cell_identify/fig1_vlnplot-new.pdf",height = 14,width = 14)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G22710')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G25710')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G79580')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G49270')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT5G57620')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G63090')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G13620')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G20840')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
plot1<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G37090')
plot2<-VlnPlot(object = wang,pt.size=0.1, features ='AT2G42840')
plot3<-VlnPlot(object = wang,pt.size=0.1, features ='AT1G77690')
plot4<-VlnPlot(object = wang,pt.size=0.1, features ='AT3G54890')
CombinePlots(plots = list(plot1,plot2,plot3,plot4),legend = 'none',ncol=1)
dev.off()