#----------------------------------------------------------------------------------------------------------
#### B cell_Trajectory inference ####
#----------------------------------------------------------------------------------------------------------

###### Monocle 3 #######
library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
### read data ###
scRNA <- readRDS("scRNA_B_Singlet.rds");
data <- GetAssayData(scRNA, assay ='RNA', slot = 'counts');
cell_metadata <- scRNA@meta.data;
gene_annotation <-data.frame(gene_short_name = rownames(data));
rownames(gene_annotation) <-rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata =cell_metadata, 
                         gene_metadata =gene_annotation);
cds <- preprocess_cds(cds, num_dim = 50);
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

### Fig.6A ###
pdf(file="monocle3_Bcells_trajectory-seuratclusterID.pdf", width=3.8, height=3.6)
plot_cells(cds,color_cells_by = "seurat_clusters",label_groups_by_cluster = TRUE, label_leaves = FALSE, 
           label_branch_points = FALSE, group_label_size=3)
while (!is.null(dev.list()))  
  dev.off()


### Fig.6B ###
cds <- order_cells(cds)
pdf(file="monocle3_Bcells_pseudotime.pdf", width=4.9, height=3.6)
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, 
           label_groups_by_cluster = TRUE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
while (!is.null(dev.list()))  
  dev.off()




### Subset cells by branch ###
### branch 3 ###
cds_branch3 <- choose_graph_segments(cds) #Fig.6D right#
cds_branch3 <- preprocess_cds(cds_branch3, num_dim = 50);
cds_branch3 <- reduce_dimension(cds_branch3, preprocess_method = "PCA")
cds_branch3.embed <- cds_branch3@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds_branch3.embed),]
cds_branch3@int_colData$reducedDims$UMAP <- int.embed
cds_branch3 <- cluster_cells(cds_branch3)
cds_branch3 <- learn_graph(cds_branch3)
cds_branch3 <- order_cells(cds_branch3)
### Screening differential genes ###
modulated_genes_branch3 <- graph_test(cds_branch3, neighbor_graph = "principal_graph", cores = 4)
write.csv(modulated_genes_branch3,"modulated_genes_branch3.csv",quote=FALSE)
genes_branch3<- row.names(subset(modulated_genes_branch3, q_value == 0 & morans_I > 0.2))

### Fig.6K (pseudotime - interested genes) ###
marker_genes_branch3 <- row.names(subset(fData(cds_branch3), 
                                  gene_short_name %in% c("CD27", "IGHD")))
pdf(file="branch3_pseudotime_genes.pdf", width=12, height=3)
plot_genes_in_pseudotime(cds_branch3[marker_genes_branch3,], color_cells_by="seurat_clusters", 
                         min_expr=0.5, ncol = 2, cell_size = 1)
while (!is.null(dev.list()))  
  dev.off()



### branch 2 ###
cds_branch2 <- choose_graph_segments(cds) #Fig.6D middle#
cds_branch2 <- preprocess_cds(cds_branch2, num_dim = 50);
cds_branch2 <- reduce_dimension(cds_branch2, preprocess_method = "PCA")
cds_branch2.embed <- cds_branch2@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds_branch2.embed),]
cds_branch2@int_colData$reducedDims$UMAP <- int.embed
cds_branch2 <- cluster_cells(cds_branch2)
cds_branch2 <- learn_graph(cds_branch2)
cds_branch2 <- order_cells(cds_branch2)
### Screening differential genes ###
modulated_genes_branch2 <- graph_test(cds_branch2, neighbor_graph = "principal_graph", cores = 4)
write.csv(modulated_genes_branch2,"modulated_genes_branch2.csv",quote=FALSE)
genes_branch2<- row.names(subset(modulated_genes_branch2, q_value == 0 & morans_I > 0.2))

### Fig.6E (pseudotime - interested genes) ###
marker_genes_branch2 <- row.names(subset(fData(cds_branch2), 
                                      gene_short_name %in% c("CD27", "IGHD")))
pdf(file="branch2_pseudotime_genes.pdf", width=12, height=3)
plot_genes_in_pseudotime(cds_branch2[marker_genes_branch2,], color_cells_by="seurat_clusters", 
                         min_expr=0.5, ncol = 2, cell_size = 1)
while (!is.null(dev.list()))  
  dev.off()



#### branch 1 ###
cds_branch1 <- choose_graph_segments(cds) #Fig.6D left#
cds_branch1 <- preprocess_cds(cds_branch1, num_dim = 50);
cds_branch1 <- reduce_dimension(cds_branch1, preprocess_method = "PCA")
cds_branch1.embed <- cds_branch1@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds_branch1.embed),]
cds_branch1@int_colData$reducedDims$UMAP <- int.embed
cds_branch1 <- cluster_cells(cds_branch1)
cds_branch1 <- learn_graph(cds_branch1)
cds_branch1 <- order_cells(cds_branch1)
### Screening differential genes ###
modulated_genes_branch1 <- graph_test(cds_branch1, neighbor_graph = "principal_graph", cores = 4)
write.csv(modulated_genes_branch1,"modulated_genes_branch1.csv",quote=FALSE)
genes_branch1<- row.names(subset(modulated_genes_branch1, q_value == 0 & morans_I > 0.2))

### Fig.6H (pseudotime - interested genes) ###
marker_genes_branch1 <- row.names(subset(fData(cds_branch1), 
                                      gene_short_name %in% c("CD27", "IGHD")))
pdf(file="branch1_pseudotime_genes.pdf", width=12, height=3)
plot_genes_in_pseudotime(cds_branch1[marker_genes_branch1,], color_cells_by="seurat_clusters", 
                         min_expr=0.5, ncol = 2, cell_size = 1)
while (!is.null(dev.list()))  
  dev.off()




### pseudotime heatmap ###
library(ClusterGVis)
### Fig.6L (branch3) ###
mat_branch3 <- pre_pseudotime_matrix(cds_obj = cds_branch3,
                                  gene_list = genes_branch3)
### kmeans ###
ck_branch3 <- clusterData(exp = mat_branch3,
                       cluster.method = "kmeans",
                       cluster.num = 3)
markGenes_branch3 <- row.names(subset(fData(cds_branch3), 
                                   gene_short_name %in% c("RPS17","CCL22","CD52","RPL9","RPL10",
                                                          "TXN","FTL","EEF1G","EEF1A1","MIF",
                                                          "PKM","ENO1","FABP5","S100A6","PRDX1",
                                                          "CD27","NME1","NHP2")))
mat_branch3_sel <- pre_pseudotime_matrix(cds_obj = cds_branch3,
                                      gene_list = markGenes_branch3)
pdf('monocle3_branch3_heatmap.pdf',height = 6,width = 6,onefile = F)
visCluster(object = ck_branch3,
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = sample(rownames(mat_branch3_sel),123,replace = T))
dev.off()



### Fig.6F (branch2) ###
mat_branch2 <- pre_pseudotime_matrix(cds_obj = cds_branch2,
                                  gene_list = genes_branch2)
### kmeans ###
ck_branch2 <- clusterData(exp = mat_branch2,
                       cluster.method = "kmeans",
                       cluster.num = 3)
markGenes_branch2 <- row.names(subset(fData(cds_branch2), 
                                   gene_short_name %in% c("HLA-A","HLA-DRB1","CD58","CD44","CD99",
                                                          "CD63","S100A10","ANXA2","LMNA","TXN",
                                                          "RPS28","FTL","HSPA6","LGALS1","VIM",
                                                          "B2M","MIR155HG","KLF6","EMP3","DNAJB1","HSPB1")))
mat_branch2_sel <- pre_pseudotime_matrix(cds_obj = cds_branch2,
                                      gene_list = markGenes_branch2)
pdf('monocle3_branch2_heatmap.pdf',height = 6,width = 6,onefile = F)
visCluster(object = ck_branch2,
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = sample(rownames(mat_branch2_sel),69,replace = T))
dev.off()



### Fig.6I (branch1) ###
mat_branch1 <- pre_pseudotime_matrix(cds_obj = cds_branch1,
                                  gene_list = genes_branch1)
### kmeans ###
ck_branch1 <- clusterData(exp = mat_branch1,
                       cluster.method = "kmeans",
                       cluster.num = 3)
markGenes_branch1 <- row.names(subset(fData(cds_branch1), 
                                   gene_short_name %in% c("B2M","CD27","CD52","CXCL1","EEF1A1",
                                                          "IGHA2","LGALS1","MMP1","RPL39","CRIP1",
                                                          "RPS3A","RPS4Y1","TMSB4X","TXNIP","RPL13")))
mat_branch1_sel <- pre_pseudotime_matrix(cds_obj = cds_branch1,
                                      gene_list = markGenes_branch1)
pdf('monocle3_branch1_heatmap.pdf',height = 6,width = 6,onefile = F)
visCluster(object = ck_branch1,
           plot.type = "heatmap",
           add.sampleanno = F,
           markGenes = sample(rownames(mat_branch1_sel),54,replace = T))
dev.off()
######### Monocel 3 end #########








###### Slingshot #######
rm(list=ls())
library(slingshot)
library(Seurat)
library(RColorBrewer)
library(grDevices)
### read data  ###
scRNA <- readRDS("scRNA_B_Singlet.rds") 
scRNA <- subset(scRNA, seurat_clusters %in% c("0","1","2","3","4","5","6","7","8","9","11","12"))
scRNA.sce <- as.SingleCellExperiment(scRNA);
class(scRNA.sce);
sim <- slingshot(scRNA.sce,clusterLabels="seurat_clusters",reducedDim="PCA");
summary(sim$slingPseudotime_1);

### Fig.6C ###
lin1 <- getLineages(sim, 
                    clusterLabels = "seurat_clusters", 
                    reducedDim = "UMAP")
pdf(file="slingshot_track-Bcells.pdf", width=9, height=9)
plot(reducedDims(sim)$UMAP,col = colorRampPalette(c("#F8766D","#0CB702","#00B8E7","#C77CFF","#FF68A1"))(13)[scRNA@meta.data[["RNA_snn_res.0.4"]]],pch=16,asp=1.3)+
  lines(SlingshotDataSet(lin1), lwd=4,col = 'black',type = 'lineages')
while (!is.null(dev.list()))  
  dev.off()



