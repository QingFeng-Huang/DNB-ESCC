#----------------------------------------------------------------------------------------------------------
#### ESCC _ dimension reduction and rename clusters ####
#----------------------------------------------------------------------------------------------------------
library(Seurat)
library(clustree)
### read data ###
data <- readRDS("scRNA_ESCC.rds");
data <- subset(x = data, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 15)
data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(object = data, 
                             selection.method = "vst", 
                             nfeatures = 2000);
data <- ScaleData(object = data, vars.to.regress = c("nCount_RNA", "percent.mt"));
data <- RunPCA(object = data,  npcs = 30, verbose = TRUE)

### Determine optimum PCA ###
sce=data
# Determine percent of variation associated with each PC   
pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Elbow plot to visualize 
ElbowPlot(data,30) 
### Overall, choose 30PC ###

### Clustree determine optimum resolution ###
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = seq(from = 0.1,to = 3.0,by = 0.1))
data <- RunUMAP(data, reduction = "pca", dims = 1:30)
p<-clustree(data, prefix = "RNA_snn_res.", node_label = "RNA_snn_res.")
### Overall, choose 2.6 resolution ###

### UMAP ###
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 2.6, algorithm = 1) 
data <- RunUMAP(object = data, dims = 1:30)
#### rename ESCC cluster ####
current.cluster.ids <-  c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
                          "21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39",
                          "40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55");
new.cluster.ids <- 
  c("T/NK cells",	"T/NK cells",	"T/NK cells",	"Fibroblasts",	"T/NK cells",	"B cells",	"Fibroblasts",
    "Plasma cells",	"T/NK cells",	"Fibroblasts",	"Myeloid cells",	"B cells",	"Myeloid cells",
    "Myeloid cells",	"Fibroblasts",	"B cells",	"Plasma cells",	"Epithelial cell",	"B cells",
    "Myeloid cells",	"Smooth muscle cells",	"Endothelial cells",	"Endothelial cells",
    "Epithelial cell",	"Epithelial cell",	"Fibroblasts",	"Mast cells",	"Myeloid cells",
    "T/NK cells",	"T/NK cells",	"T/NK cells",	"Fibroblasts",	"Cycling cells",	"Fibroblasts",
    "Epithelial cell",	"Endothelial cells",	"pDC",	"Epithelial cell",	"Fibroblasts",
    "T/NK cells",	"B cells",	"Fibroblasts",	"B cells",	"Myeloid cells",	"Epithelial cell",
    "Cycling cells",	"B cells",	"Epithelial cell",	"B cells",	"Epithelial cell",
    "Epithelial cell",	"Mast cells",	"Plasma cells",	"Fibroblasts",	"Plasma cells",	"Endothelial cells");
names(new.cluster.ids) <- levels(data);
data <- RenameIdents(data, new.cluster.ids);
data$celltype <- data@active.ident
saveRDS(data, "scRNA_ESCC_newcluster.rds")



####### extract B cell, T&NK cell ########
part = scRNA[, Idents(scRNA) %in% c("5","11","15","18","40","46","48","42")]
saveRDS(part, file = "scRNA_B.rds");
part = scRNA[, Idents(scRNA) %in% c("0","1","2","4","8","28","29","30","39")]
saveRDS(part, file = "scRNA_T&NK.rds")




#----------------------------------------------------------------------------------------------------------
#### Fig.S5A,S5B,S5C ####
#----------------------------------------------------------------------------------------------------------
rm(list=ls())
library(Seurat)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)
###### read data ######
data <- readRDS("scRNA_ESCC_newcluster.rds");

### Fig.S5A (UMAP) ###
p1<- DimPlot(data, raster=FALSE, reduction = "umap", group.by = "ident", 
             pt.size=0.2,label = TRUE, label.size = 5.5, repel = TRUE)+
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank())
ggsave(filename = "UMAP_ESCC_newcluster.pdf", plot = p1, device = 'pdf', width = 9, height = 6, units = 'in')


### Fig.S5B (UMAP - interested genes) ###
pdf(file="ESCC Marker.pdf",width=21, height=20)
FeaturePlot(data, features = c("KRT17","KRT8","CLDN5","DCN",
                               "TAGLN","LILRA4","LYZ","CD68",
                               "TPSAB1","CD3D","TRBC1","CD79A",
                               "MS4A1","MZB1","MKI67"))
while (!is.null(dev.list()))     dev.off()


### Fig.S5C (z-score Heatmap) ###
Idents(data)= data@meta.data[["seurat_clusters"]]
data2 <- GetAssayData(data,slot="counts",assay="RNA");
data2 <- GetAssayData(data,slot="data",assay="RNA");
data2 <- GetAssayData(data,slot="scale.data",assay="RNA");
gene <- read.table('genes.txt',header=T,sep="\t",quote='',stringsAsFactors = F)
marker_genes <- AverageExpression(data,
                                  assays = "RNA",
                                  features = gene$gene,
                                  group.by = "celltype",
                                  slot="data");
marker_genes_info <- marker_genes$RNA;
write.table(marker_genes_info,"ESCC_marker_genes_info.csv",quote=F,row.names=TRUE,sep=",")
###### after sorting items manually, read table ######
data3 <- read.csv("ESCC_marker_genes_info_order.csv",header = T,row.names = 1,check.names = F)
exp <- apply(data3, 1, scale)
rownames(exp) <- colnames(data3);
exp <- t(exp);
col_fun =colorRamp2(c(-2, 0 ,2), c("#3333FF","white","#FF0000"))
pdf(file="ESCC_markers_z-score_heatmap.pdf", 
    width=4.5, height=7)
Heatmap(exp, name = "z-score", 
        column_order = colnames(exp),
        row_order = rownames(exp),
        col = col_fun)
while (!is.null(dev.list()))  
  dev.off()

