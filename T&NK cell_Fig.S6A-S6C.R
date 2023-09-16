#----------------------------------------------------------------------------------------------------------
#### T&NK cell _ dimension reduction and find double cell ####
#----------------------------------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
### read data ###
data <- readRDS("scRNA_T&NK.rds");
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
ElbowPlot(data,30)#综上选择20PC
### Overall, choose 20PC ####
### UMAP ####
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data, resolution = 0.9, algorithm = 1) 
data <- RunUMAP(object = data, dims = 1:20)
### find double cell ###
scRNA <- data 
sweep.res.list <- paramSweep_v3(scRNA, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() 
DoubletRate = ncol(scRNA)*8*1e-6 
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) 
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scRNA <- doubletFinder_v3(scRNA, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
head(x = scRNA[[]]);
table(scRNA$DF.classifications_0.25_0.25_5980);
scRNA2 <- subset(scRNA,DF.classifications_0.25_0.25_5980 == "Singlet")
saveRDS(scRNA2, file = "scRNA_T&NK_Singlet.rds")



#----------------------------------------------------------------------------------------------------------
#### Fig.S6A,S6B,S6C ####
#----------------------------------------------------------------------------------------------------------
rm(list=ls())
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
### read data  ###
data <- readRDS("scRNA_T&NK_Singlet.rds");

### Fig.S6A (UMAP) ###
p1 <- DimPlot(object = data, reduction = "umap", label = TRUE);
ggsave(filename = "UMAP_T&NK cell.pdf", plot = p1, width = 6, height = 5)


### Fig.S6C (UMAP - interested genes) ###
pdf(file="T&NK cell Marker.pdf", 
    width=18, height=12)
FeaturePlot(data, features = c("CD3E","CD4","CD8A",
                               "CCR7","FOXP3","TYROBP",
                               "CTLA4","CD69","NKG7",
                               "KLRC1","CMC1","ANKRD28"))
while (!is.null(dev.list()))     dev.off()


### Fig.S6B (z-score Heatmap) ###
Idents(data)= data@meta.data[["seurat_clusters"]]
data2 <- GetAssayData(data,slot="counts",assay="RNA");
data2 <- GetAssayData(data,slot="data",assay="RNA");
data2 <- GetAssayData(data,slot="scale.data",assay="RNA");
gene <- read.table('T&NKcellgenes.txt',header=T,sep="\t",quote='',stringsAsFactors = F)
marker_genes <- AverageExpression(data,
                                  assays = "RNA",
                                  features = gene$gene,
                                  group.by = "seurat_clusters",
                                  slot="data");
marker_genes_info <- marker_genes$RNA;
write.table(marker_genes_info,"T&NKcell_marker_genes_info.csv",quote=F,row.names=TRUE,sep=",")
### after sorting items manually, read table ###
data3 <- read.csv("T&NKcell_marker_genes_info_order.csv",header = T,row.names = 1,check.names = F)
exp <- apply(data3, 1, scale)
rownames(exp) <- colnames(data3);
exp <- t(exp);
col_fun =colorRamp2(c(-2, 0 ,2), c("#3333FF","white","#FF0000"))
pdf(file="T&NKcell_markers_z-score_heatmap.pdf", 
    width=4.5, height=5)
Heatmap(exp, name = "z-score", 
        column_order = colnames(exp),
        row_order = rownames(exp),
        col = col_fun)
while (!is.null(dev.list()))  
  dev.off()
######## end #########



#----------------------------------------------------------------------------------------------------------
#### rename T&NK cell cluster ####
#----------------------------------------------------------------------------------------------------------
current.cluster.ids <-  c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19");
new.cluster.ids <- 
  c("CD4 Tn",	"CMC1+ CD8 T",	"Treg",	"CD8 T",	"CD69+ CD8 Tmem",	"CMC1+ CD8 T",	"CD4 Tem","NKcyto", 
    "KLRC1+ CD8 T",	"Treg",	"CD69+ CD8 Tmem",	"CD4 Tem",	"CD8 Tex",	"CD69+ CD8 Tmem", "CD4 Tex", "CD69+ CD8 Tmem",
    "NKrest", "CD8 T", "CD8 T", "CD4 Tex");
names(new.cluster.ids) <- levels(data);
T_scRNA <- RenameIdents(data, new.cluster.ids);
T_scRNA$celltype <- T_scRNA@active.ident
saveRDS(T_scRNA, "T&NK cell_newcluster.rds")



