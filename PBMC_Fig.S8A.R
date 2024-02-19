#----------------------------------------------------------------------------------------------------------
#### PBMC _ dimension reduction ####
#----------------------------------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
### read data ###
data <- readRDS("scRNA_PBMC.rds");
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
### Overall, choose 20PC ####
### UMAP ####
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data, resolution = 0.1, algorithm = 1) 
data <- RunUMAP(object = data, dims = 1:20)



#----------------------------------------------------------------------------------------------------------
#### Fig.S8A  ####
#----------------------------------------------------------------------------------------------------------
### Fig.S8A (UMAP) ###
p1 <- DimPlot(object = data, reduction = "umap", label = TRUE);
ggsave(filename = "UMAP_PBMC.pdf", plot = p1, width = 6, height = 5)


### Fig.S8A (UMAP - interested genes) ###
pdf(file="PBMC Marker.pdf", 
    width=13.5, height=12)
FeaturePlot(data, features = c("CD79A","MS4A1","MZB1","CD3D",
                               "TYROBP","KLRB1","LYZ","LILRA4"))
while (!is.null(dev.list()))     dev.off()



####### extract B cell & Plasma cell #######
part = data[, Idents(data) %in% c("2","5")]
saveRDS(part, file = "scRNA_PBMC_B&Plasma.rds")


