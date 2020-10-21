library(Seurat)
library(dplyr)
 library(magrittr)
 library(gtools)
 library(stringr)
 library(Matrix)
 setwd("D://data/ScRNAcode")
##读入STAR数据
 matrix.dir="STAR/"
 barcode.path <- paste0(matrix.dir,"barcodes.tsv")
 features.path <- paste0(matrix.dir,"features.tsv")
 matrix.path <- paste0(matrix.dir, "matrix.mtx")
 STARmatrix <- readMM(file = matrix.path)
 feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
 barcode.names = read.delim(barcode.path,
                              header = FALSE,
                             stringsAsFactors = FALSE)
 colnames(STARmatrix) = barcode.names$V1
 rownames(STARmatrix) = feature.names$V2
 STARmatrix<-as.matrix(STARmatrix)
 STARmatrix[1:6,1:6]
##创建seurat对象
##创建STAR的Seurat对象
 STAR <- CreateSeuratObject(STARmatrix,project = "zsz",
                             min.cells = 3, min.features = 200)
##创建Cellranger的Seurat对象
 dir="cellranger/"
 counts <- Read10X(data.dir = dir)
 RANGER <- CreateSeuratObject(counts, project = "zsz", 
                               min.cells=3, min.features = 200)
 
##数据比较
dim(STARmatrix)
dim(counts)
dim(STAR)
dim(RANGER)
fivenum(apply(STARmatrix,1,function(x) sum(x0)))
fivenum(apply(counts,1,function(x) sum(x0)))
pdf("box.pdf",height = 9,width = 9)
boxplot(apply(STARmatrix,1,function(x) sum(x0) ),main = "STAR",col = "lightgray")
boxplot(apply(counts,1,function(x) sum(x0) ),main = "Cellranger",col = "lightgray")
dev.off()
pdf("hist.pdf",height = 9,width = 9)
hist(apply(STARmatrix,2,function(x) sum(x0) ),col = "lightgray",
     breaks=20,xlim=c(0,4000),ylim=c(0,800),
     labels=F,main="STAR",
     xlab="genes",ylab="cells")
abline(v=median(apply(STARmatrix,2,function(x) sum(x0))),col='red')
hist(apply(counts,2,function(x) sum(x0) ),col = "lightgray",
     breaks=20,xlim=c(0,4000),ylim=c(0,800),
     labels=F,main="Cellranger",
     xlab="genes",ylab="cells")
abline(v=median(apply(counts,2,function(x) sum(x0))),col='red')
dev.off()              

pdf("qc.pdf",height = 9,width = 9)
VlnPlot(STAR,
        features = c("nFeature_RNA", "nCount_RNA"), 
        pt.size = 0.1,
        ncol = 2)
VlnPlot(RANGER,
        features = c("nFeature_RNA", "nCount_RNA"), 
        pt.size = 0.1, 
        ncol = 2)
dev.off()

STAR<-subset(STAR,subset=nFeature_RNA>500 & nFeature_RNA<2000)
STAR <- NormalizeData(STAR, normalization.method = "LogNormalize", scale.factor = 10000)
table(STAR@meta.data$orig.ident)
STAR <- FindVariableFeatures(STAR, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(STAR), 10) 
top10

RANGER<-subset(RANGER,subset=nFeature_RNA>500 & nFeature_RNA<2000)
RANGER<- NormalizeData(RANGER, normalization.method = "LogNormalize", scale.factor = 10000)
table(RANGER@meta.data$orig.ident)
RANGER <- FindVariableFeatures(RANGER, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RANGER), 10) 
top10

scale.genes <-  rownames(STAR)
STAR <- ScaleData(STAR, features = scale.genes)
STAR <- RunPCA(STAR, features = VariableFeatures(STAR)) 
plot1 <- ElbowPlot(STAR, ndims=30, reduction="pca") 
scale.genes <-  rownames(RANGER)
RANGER <- ScaleData(RANGER, features = scale.genes)
RANGER <- RunPCA(RANGER, features = VariableFeatures(RANGER)) 
plot2 <- ElbowPlot(RANGER, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("pca.pdf", plot = plotc, width = 8, height = 4)


STAR <- FindNeighbors(STAR, dims = 1:10) 
STAR <- FindClusters(STAR, resolution = 0.8)
table(STAR@meta.data$seurat_clusters)
metadata <- STAR@meta.data
cell_cluster <-data.frame(cell_ID=rownames(metadata),
                          cluster_ID=metadata$seurat_clusters)
STAR <- RunUMAP(STAR, dims = 1:20)
embed_tsne <- Embeddings(STAR, 'umap')
plot1 = DimPlot(STAR, reduction = "umap" ,label = "T", pt.size = 1,label.size = 4)
RANGER <- FindNeighbors(RANGER, dims = 1:10) 
RANGER <- FindClusters(RANGER, resolution = 0.8)
table(RANGER@meta.data$seurat_clusters)
metadata <- RANGER@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
RANGER <- RunUMAP(RANGER,n.neighbors = 30,dims = 1:20)
embed_umap <- Embeddings(RANGER, 'umap')
plot2 = DimPlot(RANGER, reduction = "umap" ,label = "T", pt.size = 1,label.size = 4)
plotc <- plot1+plot2
ggsave("umap.pdf", plot = plotc, width = 8, height = 4)

