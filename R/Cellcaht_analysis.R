library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(CellChat)
library(ggplot2)
library(ggalluvial)
obj.combined <-readRDS("./Analysis/data_resolu0.6.rds")
current.cluster.ids <- c("RB1-rep1","RB1-rep2","RB2-rep1","RB2-rep2","RB3-rep1","RB3-rep2","RB4","RB5","RB6","RB7")
new.cluster.ids <- c("RB1","RB1","RB2","RB2","RB3","RB3","RB4","RB5","RB6","RB7")
obj.combined@meta.data$group <- plyr::mapvalues(x = obj.combined@meta.data[,"sample"], from = current.cluster.ids, to = new.cluster.ids)
Mac_Cluster15_obj.combined <-subset(x = obj.combined,idents="15")
Mac_Cluster15_obj.combined
DefaultAssay(Mac_Cluster15_obj.combined) <- "integrated"
#Cluster the cells 
Mac_Cluster15_obj.combined <- FindNeighbors(Mac_Cluster15_obj.combined,reduction = "pca",dims = 1:10)
Mac_Cluster15_obj.combined <- FindClusters(Mac_Cluster15_obj.combined,resolution = 0.3)
Mac_Cluster15_obj.combined <- RunTSNE(Mac_Cluster15_obj.combined, reduction = "pca", dims = 1:10)
Mac_Cluster15_obj.combined <- RunUMAP(Mac_Cluster15_obj.combined, reduction = "pca", dims = 1:10)

bb <-as.matrix(obj.combined@meta.data)
aa <-as.matrix(Mac_Cluster15_obj.combined@meta.data)
for(i in 1:nrow(aa)){
bb[which(rownames(bb)==rownames(aa)[i]),13] <-paste0('Glia_',aa[i,17])
}
obj.combined@meta.data <-data.frame(bb)

sub_obj.combined <-subset(obj.combined,idents=c("0","1","2","3","5","6","7","8","9","10","11","12","13","Glia_0","Glia_1","16"))
head(x = obj.combined@meta.data)
current.cluster.ids <- c("0","1","2","5","6","3","7","8","9","10","13","16","11","12","Glia_1","Glia_0")
new.cluster.ids <- c("CP-C0","CP-C1","CP-C2","CP-C5","CP-C6","MKI67+ CP-C3","MKI67+ CP-C7","MKI67+ CP-C8","MKI67+ CP-C9","MKI67+ CP-C10","Rod-like-C13","Cone-like-C16","Neural cell-C11","CAF-C12","Astrocyte","Macrophage")
sub_obj.combined@meta.data$celltype <- plyr::mapvalues(x = sub_obj.combined@meta.data[,13], from = current.cluster.ids, to = new.cluster.ids)
head(x = sub_obj.combined@meta.data)
sub_obj.combined@active.ident <-factor(as.matrix(sub_obj.combined@meta.data)[,17])
Idents(sub_obj.combined) <- factor(Idents(sub_obj.combined), levels = c("CP-C0","CP-C1","CP-C2","CP-C5","CP-C6","MKI67+ CP-C3","MKI67+ CP-C7","MKI67+ CP-C8","MKI67+ CP-C9","MKI67+ CP-C10","Cone-like-C16","Rod-like-C13","Neural cell-C11","CAF-C12","Astrocyte","Macrophage"))

setwd("./cellchat")
cellchat <- createCellChat(sub_obj.combined@assays$RNA@data) ##creat cellchat object
identity = data.frame(group =sub_obj.combined@active.ident, row.names = names(sub_obj.combined@active.ident)) #  cell labels
unique(identity$group) # check the cell labels
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels") ##add lables
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human  ##input human database
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object
##CellChatDB$interaction[1:4,1:4]
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 8) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use=TRUE)
cellchat <- filterCommunication(cellchat,min.cells=10)
#cellchat <- filterCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
###
#cellchat@netP$pathways
#head(cellchat@LR$LRsig)
##可视化
pathways.show <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
pathways.show
##
for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show[i], layout = "circle", vertex.weight = groupSize,pt.title=20,vertex.label.cex = 1.7)
#  netVisual(cellchat, signaling = pathways.show[i], layout = "chord", vertex.weight = groupSize,pt.title=20,vertex.label.cex = 1.7)
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  # Visualize signaling roles of cell groups
  pdf(file = paste0(pathways.show[i], "_signalRole.pdf"))
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.off()
}
##
for (i in 1:length(pathways.show)) {
pdf(file =paste0(pathways.show[i], "_chord_adjust.pdf"))
group.cellType <- c("CP-C0","CP-C1","CP-C2","CP-C5","CP-C6","MKI67+ CP-C3","MKI67+ CP-C7","MKI67+ CP-C8","MKI67+ CP-C9","MKI67+ CP-C10","Rod-like-C13","Cone-like-C16","Neural cell-C11","CAF-C12","Astrocyte","Macrophage") 
names(group.cellType) <- levels(cellchat@idents) 
netVisual_chord_cell(cellchat, signaling = pathways.show[i], group = group.cellType, title.name = paste0(pathways.show[i], " signaling network")) 
dev.off()
}
for (i in 1:length(pathways.show)) {
  pdf(file = paste0(pathways.show[i], "_heatmap.pdf"))
  netVisual_heatmap(cellchat, signaling = pathways.show[i], color.heatmap = "Reds") 
  dev.off()
} 
##
setwd("./cellchat/function")
ident <-c(1:16)
for (i in 1:length(ident)) {
pdf(file =paste0(levels(cellchat@idents)[i], "_LR.pdf"))
print(netVisual_bubble(cellchat, sources.use = i, targets.use = ident[-i], remove.isolate = FALSE))
dev.off()
}
ident <-c(1:16)
for (i in 1:length(ident)) {
pdf(file =paste0(levels(cellchat@idents)[i], "_netVisual_chord.pdf"))
netVisual_chord_gene(cellchat, sources.use = i, targets.use = ident[-i], lab.cex = 0.5,legend.pos.y = 30) 
dev.off()
}
# show all the significant interactions (L-R pairs) associated with certain signaling pathways 
pdf(file ="test_netVisual_chord.pdf")
netVisual_chord_gene(cellchat, sources.use = c(14:16), targets.use = c(1:13),legend.pos.x = 8) 
dev.off()
#> Note: The second link end is drawn out of sector 'CXCR4 '. #> Note: The first link end is drawn out of sector 'CXCL12 '.
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use') 
pdf(file ="test_netVisual_chord.pdf")
netVisual_chord_gene(cellchat, sources.use = c(14:16), targets.use = c(1:13), slot.name = "netP", legend.pos.x = 10) 
dev.off()
#> Note: The second link end is drawn out of sector ' '. #> Note: The first link end is drawn out of sector 'MIF'. #> Note: The second link end is drawn out of sector ' '. #> Note: The first link end is drawn out of sector 'CXCL '.
###
pdf(file = "Communication_strength_weight.pdf")
par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions") 
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
pdf(file = "Communication_weight.pdf")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE) 
for (i in 1:nrow(mat)) { 
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)) 
mat2[i, ] <- mat[i, ] 
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]) 
}
dev.off()
##
nPatterns = 5
pdf(file = paste0("CommunicationPatterns_sender_heatmap.pdf"))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()
pdf(file = "patternAnalysis_sender_river.pdf", width = 7, height = 4)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
gg <- netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave(filename="patternAnalysis_sender_dot.pdf", plot=gg, width = 5.5, height = 4, units = 'in', dpi = 300)
pdf(file = paste0("CommunicationPatterns_receiver_heatmap.pdf"))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()
pdf(file = "patternAnalysis_receiver_river.pdf", width = 7, height = 4)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
gg <- netAnalysis_dot(cellchat, pattern = "incoming")
ggsave(filename="patternAnalysis_receiver_dot.pdf", plot=gg, width = 5.5, height = 4, units = 'in', dpi = 300)
##
cellchat <- computeNetSimilarity(cellchat, type = "functional", thresh = 0.25)    
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional", k = 4)
gg <- netVisual_embedding(cellchat, type = "functional", pathway.remove.show = F)
cowplot::save_plot("2Dmanifold_FunctionalSimilarity_signalingPathways.pdf", gg, base_height = 3, base_width = 4)
pdf(file = "2Dmanifold_FunctionalSimilarity_signalingPathways_zoomIn.pdf", width = 2, height = 2.5*3)
netVisual_embeddingZoomIn(cellchat, type = "functional")
dev.off()
cellchat <- computeNetSimilarity(cellchat, type = "structural", thresh = 0.25)
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
gg <- netVisual_embedding(cellchat, type = "structural")
cowplot::save_plot("2Dmanifold_StructureSimilarity_signalingPathways.pdf", gg, base_height = 3, base_width = 4)
pdf(file = "2Dmanifold_StructureSimilarity_signalingPathways_zoomIn.pdf", width = 2, height = 2.5*3)
netVisual_embeddingZoomIn(cellchat, type = "structural")
dev.off()
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing") 
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming") 
pdf(file = "heatmap_outgoing_incoming.pdf", width =10, height = 5)
print(ht1 + ht2)
dev.off()
##gene
setwd("./cellchat/function/gene")
for (i in 1:length(pathways.show)) 
{
pdf(file = paste0(pathways.show[i], "_gene.pdf"))
print(plotGeneExpression(cellchat, signaling =pathways.show[i]))
dev.off()
}
saveRDS(cellchat,"Cellchat.rds")

