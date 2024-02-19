#----------------------------------------------------------------------------------------------------------
#### DN B, Memory B, Naive B and T&NK cell_Cell-cell communication ####
#----------------------------------------------------------------------------------------------------------
### Fig.4G-3I, Supplementary Table S6 ###
### CellChat ###
library(CellChat)
### read data ###
scRNA_T <- readRDS("scRNA_B_Singlet.rds");
scRNA_B <- readRDS("scRNA_T&NK_Singlet.rds");
unique(scRNA_B@meta.data[["celltype"]])
scRNA_B <- subset(scRNA_B,celltype %in% c ("DN B","Memory B","Naive B"));
mergeCols <- c("RNA_snn_res.2.6")
scRNA <- merge(scRNA_T, scRNA_B, by = mergeCols, all = TRUE);
scRNA;
head(x=scRNA[[]]);
View(scRNA);
unique(scRNA@meta.data[["celltype"]]);
unique(scRNA@meta.data[["seurat_clusters"]]);
scRNA@commands$FindClusters
data.input  <- scRNA@assays$RNA@data
identity = data.frame(group = scRNA@meta.data[["celltype"]], 
                      row.names = names(scRNA@meta.data[["celltype"]])) 
unique(identity$group)
cellchat <- createCellChat(object = data.input)
cellchat
summary(cellchat)
str(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling"))
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
options(future.globals.maxSize = 1000*1024^2)
cellchat <- computeCommunProb(cellchat)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"cellchat_net_lr.csv",quote=F,row.names=FALSE) # Supplementary Table S3 #
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name="netP")
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
levels(cellchat@idents) 
vertex.receiver = seq(1,4)
cellchat@LR$LRsig$pathway_name
cellchat@LR$LRsig$antagonist
pathways_name <- cellchat@netP$pathways;
pdf(file=paste("","Number of interactions.pdf"), 
    width=27, height=12)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale=T,label.edge=F,
                 arrow.width = 1,
                 arrow.size = 0.7,
                 title.name="Number of interactions")
while (!is.null(dev.list()))  
  dev.off()
mat <- cellchat@net$count
allname <- rownames(mat)
par(mfow = c(3,4),xpd = TRUE)
for (j in 1:nrow(mat)){
  mat2 <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
  mat2[j, ]<-mat[j, ]
  name=allname[j]
  pdf(file=paste("",name,"_Number of interactions",".pdf"), 
      width=27, height=12)
  p<-netVisual_circle(mat2,vertex.weight = groupSize,weight.scale=T,edge.weight.max=max(mat),
                      arrow.width = 1,
                      arrow.size = 0.7,
                      title.name=rownames(mat)[j])
  while (!is.null(dev.list()))  
    dev.off()
  print(j)
}


