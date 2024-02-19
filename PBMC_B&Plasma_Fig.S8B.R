#----------------------------------------------------------------------------------------------------------
#### PBMC_B&Plasma _ dimension reduction and find double cell ####
#----------------------------------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
###### read data #########
data <- readRDS("scRNA_PBMC_B&Plasma.rds");
data <- subset(x = data, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 15)
data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(object = data, 
                             selection.method = "vst", 
                             nfeatures = 2000);
data <- ScaleData(object = data, vars.to.regress = c("nCount_RNA", "percent.mt"));
data <- RunPCA(object = data,  npcs = 30, verbose = TRUE)
#### Determine optimum PCA ####
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
### Overall, choose 15PC ####
### UMAP ####
data <- FindNeighbors(data, reduction = "pca", dims = 1:15)
data <- FindClusters(data, resolution = 0.2, algorithm = 1) 
data <- RunUMAP(object = data, dims = 1:15)
### find double cell ###
scRNA <- data 
sweep.res.list <- paramSweep_v3(scRNA, PCs = 1:15, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() 
DoubletRate = ncol(scRNA)*8*1e-6 
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) 
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scRNA <- doubletFinder_v3(scRNA, PCs = 1:15, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
head(x = scRNA[[]]);
table(scRNA$DF.classifications_0.25_0.08_82);
scRNA2 <- subset(scRNA,DF.classifications_0.25_0.08_82 == "Singlet")
saveRDS(scRNA2, file = "scRNA_PBMC_B&Plasma_Singlet.rds")



#----------------------------------------------------------------------------------------------------------
#### Fig.S8B ####
#----------------------------------------------------------------------------------------------------------
rm(list=ls())
library(Seurat)
library(ggplot2)
###### read data  ######
data <- readRDS("scRNA_PBMC_B&Plasma_Singlet.rds");

### Fig.S8B (UMAP) ###
p2 <- DimPlot(object = data, reduction = "umap", label = TRUE);
ggsave(filename = "UMAP_PBMC_B&Plasma.pdf", plot = p2, width = 4.5, height = 3.6)


### Fig.S8B (Violin - interested genes) ###
pdf(file="PBMC_B&Plasma Marker-Violin.pdf", 
    width=9, height=8)
VlnPlot(data, features = c("MS4A1","CD27","IGHD","MZB1"), pt.size=0, slot = "data")
while (!is.null(dev.list()))     dev.off()

