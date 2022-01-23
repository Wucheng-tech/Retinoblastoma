library(Seurat)
library(dplyr)
setwd("./Analysis")
######Setup the Seurat objects#######
RB1_1.data = Read10X("./data/RB1_rep1/RB1_rep1/outs/filtered_feature_bc_matrix")
RB1_2.data = Read10X("./data/RB1_rep2/RB1_rep2/outs/filtered_feature_bc_matrix")
RB2_1.data = Read10X("./data/RB2_rep1/RB2_rep1/outs/filtered_feature_bc_matrix")
RB2_2.data = Read10X("./data/RB2_rep2/RB2_rep2/outs/filtered_feature_bc_matrix")
RB3_1.data = Read10X("./data/RB3_rep1/RB3_rep1/outs/filtered_feature_bc_matrix")
RB3_2.data = Read10X("./data/RB3_rep2/RB3_rep2/outs/filtered_feature_bc_matrix")
RB4.data = Read10X("./data/RB4/RB4_rep1/outs/filtered_feature_bc_matrix")
RB5.data = Read10X("./data/RB5/RB5_rep1/outs/filtered_feature_bc_matrix")
RB6.data = Read10X("./data/RB6/RB6_rep1/outs/filtered_feature_bc_matrix")
RB7.data = Read10X("./data/RB7/RB7_rep1/outs/filtered_feature_bc_matrix")


obj1_1 = CreateSeuratObject(RB1_1.data, project = "RB1_1", min.cells=10)
obj1_2 = CreateSeuratObject(RB1_2.data, project = "RB1_2", min.cells=10)
obj2_1 = CreateSeuratObject(RB2_1.data, project = "RB2_1", min.cells=10)
obj2_2 = CreateSeuratObject(RB2_2.data, project = "RB2_2", min.cells=10)
obj3_1 = CreateSeuratObject(RB3_1.data, project = "RB3_1", min.cells=10)
obj3_2 = CreateSeuratObject(RB3_2.data, project = "RB3_2", min.cells=10)
obj4 = CreateSeuratObject(RB4.data, project = "RB4", min.cells=10)
obj5 = CreateSeuratObject(RB5.data, project = "RB5", min.cells=10)
obj6 = CreateSeuratObject(RB6.data, project = "RB6", min.cells=10)
obj7 = CreateSeuratObject(RB7.data, project = "RB7", min.cells=10)

obj1_1$sample = "RB1-rep1"
obj1_2$sample = "RB1-rep2"
obj2_1$sample = "RB2-rep1"
obj2_2$sample = "RB2-rep2"
obj3_1$sample = "RB3-rep1"
obj3_2$sample = "RB3-rep2"
obj4$sample = "RB4"
obj5$sample = "RB5"
obj6$sample = "RB6"
obj7$sample = "RB7"

obj1_1[["percent.mt"]]= PercentageFeatureSet(obj1_1, pattern = "^(M|m)(T|t)-")
obj1_2[["percent.mt"]]= PercentageFeatureSet(obj1_2, pattern = "^(M|m)(T|t)-")
obj2_1[["percent.mt"]]= PercentageFeatureSet(obj2_1, pattern = "^(M|m)(T|t)-")
obj2_2[["percent.mt"]]= PercentageFeatureSet(obj2_2, pattern = "^(M|m)(T|t)-")
obj3_1[["percent.mt"]]= PercentageFeatureSet(obj3_1, pattern = "^(M|m)(T|t)-")
obj3_2[["percent.mt"]]= PercentageFeatureSet(obj3_2, pattern = "^(M|m)(T|t)-")
obj4[["percent.mt"]]= PercentageFeatureSet(obj4, pattern = "^(M|m)(T|t)-")
obj5[["percent.mt"]]= PercentageFeatureSet(obj5, pattern = "^(M|m)(T|t)-")
obj6[["percent.mt"]]= PercentageFeatureSet(obj6, pattern = "^(M|m)(T|t)-")
obj7[["percent.mt"]]= PercentageFeatureSet(obj7, pattern = "^(M|m)(T|t)-")

obj1_1[["percent.RP"]]= PercentageFeatureSet(obj1_1, pattern = "^RPL|^RPS")
obj1_2[["percent.RP"]]= PercentageFeatureSet(obj1_2, pattern = "^RPL|^RPS")
obj2_1[["percent.RP"]]= PercentageFeatureSet(obj2_1, pattern = "^RPL|^RPS")
obj2_2[["percent.RP"]]= PercentageFeatureSet(obj2_2, pattern = "^RPL|^RPS")
obj3_1[["percent.RP"]]= PercentageFeatureSet(obj3_1, pattern = "^RPL|^RPS")
obj3_2[["percent.RP"]]= PercentageFeatureSet(obj3_2, pattern = "^RPL|^RPS")
obj4[["percent.RP"]]= PercentageFeatureSet(obj4, pattern = "^RPL|^RPS")
obj5[["percent.RP"]]= PercentageFeatureSet(obj5, pattern = "^RPL|^RPS")
obj6[["percent.RP"]]= PercentageFeatureSet(obj6, pattern = "^RPL|^RPS")
obj7[["percent.RP"]]= PercentageFeatureSet(obj7, pattern = "^RPL|^RPS")
#obj10[["percent.RP"]]= PercentageFeatureSet(obj10, pattern = "^(M|m)(T|t)-")


HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
HB_m <- match(HB.genes_total,rownames(obj1_1@assays$RNA))
HB.genes <- rownames(obj1_1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj1_1[["percent.HB"]]<-PercentageFeatureSet(obj1_1,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj1_2@assays$RNA))
HB.genes <- rownames(obj1_2@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj1_2[["percent.HB"]]<-PercentageFeatureSet(obj1_2,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj2_1@assays$RNA))
HB.genes <- rownames(obj2_1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj2_1[["percent.HB"]]<-PercentageFeatureSet(obj2_1,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj2_2@assays$RNA))
HB.genes <- rownames(obj2_2@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj2_2[["percent.HB"]]<-PercentageFeatureSet(obj2_2,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj3_1@assays$RNA))
HB.genes <- rownames(obj3_1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj3_1[["percent.HB"]]<-PercentageFeatureSet(obj3_1,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj3_2@assays$RNA))
HB.genes <- rownames(obj3_2@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj3_2[["percent.HB"]]<-PercentageFeatureSet(obj3_2,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj4@assays$RNA))
HB.genes <- rownames(obj4@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj4[["percent.HB"]]<-PercentageFeatureSet(obj4,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj5@assays$RNA))
HB.genes <- rownames(obj5@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj5[["percent.HB"]]<-PercentageFeatureSet(obj5,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj6@assays$RNA))
HB.genes <- rownames(obj6@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj6[["percent.HB"]]<-PercentageFeatureSet(obj6,features=HB.genes)
HB_m <- match(HB.genes_total,rownames(obj7@assays$RNA))
HB.genes <- rownames(obj7@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
obj7[["percent.HB"]]<-PercentageFeatureSet(obj7,features=HB.genes)

pdf("QC.pdf",width=15,height=8)
print(VlnPlot(object = obj1_1, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj1_2, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj2_1, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj2_2, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj3_1, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj3_2, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj4, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj5, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj6, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
print(VlnPlot(object = obj7, features = c("nFeature_RNA", "percent.mt","percent.RP","percent.HB"), ncol = 4))
dev.off()


obj1_1 = subset(obj1_1, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj1_2 = subset(obj1_2, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj2_1 = subset(obj2_1, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj2_2 = subset(obj2_2, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj3_1 = subset(obj3_1, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj3_2 = subset(obj3_2, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40)
obj4 = subset(obj4, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj5 = subset(obj5, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj6 = subset(obj6, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)
obj7 = subset(obj7, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & percent.RP < 40 & percent.HB < 10)

#Normalizing the data
obj1_1 <- NormalizeData(obj1_1)
obj1_2 <- NormalizeData(obj1_2)
obj2_1 <- NormalizeData(obj2_1)
obj2_2 <- NormalizeData(obj2_2)
obj3_1 <- NormalizeData(obj3_1)
obj3_2 <- NormalizeData(obj3_2)
obj4 <- NormalizeData(obj4)
obj5 <- NormalizeData(obj5)
obj6 <- NormalizeData(obj6)
obj7 <- NormalizeData(obj7)

#Identification of highly variable features (feature selection)
obj1_1 <- FindVariableFeatures(obj1_1, selection.method = "vst", nfeatures = 2000)
obj1_2 <- FindVariableFeatures(obj1_2, selection.method = "vst", nfeatures = 2000)
obj2_1 <- FindVariableFeatures(obj2_1, selection.method = "vst", nfeatures = 2000)
obj2_2 <- FindVariableFeatures(obj2_2, selection.method = "vst", nfeatures = 2000)
obj3_1 <- FindVariableFeatures(obj3_1, selection.method = "vst", nfeatures = 2000)
obj3_2 <- FindVariableFeatures(obj3_2, selection.method = "vst", nfeatures = 2000)
obj4 <- FindVariableFeatures(obj4, selection.method = "vst", nfeatures = 2000)
obj5 <- FindVariableFeatures(obj5, selection.method = "vst", nfeatures = 2000)
obj6 <- FindVariableFeatures(obj6, selection.method = "vst", nfeatures = 2000)
obj7 <- FindVariableFeatures(obj7, selection.method = "vst", nfeatures = 2000)
#obj10 <- FindVariableFeatures(obj10, selection.method = "vst", nfeatures = 2000)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj1_1 <- CellCycleScoring(obj1_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj1_2 <- CellCycleScoring(obj1_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj2_1 <- CellCycleScoring(obj2_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj2_2 <- CellCycleScoring(obj2_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj3_1 <- CellCycleScoring(obj3_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj3_2 <- CellCycleScoring(obj3_2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj4 <- CellCycleScoring(obj4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj5 <- CellCycleScoring(obj5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj6 <- CellCycleScoring(obj6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj7 <- CellCycleScoring(obj7, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

########Perform integration###############
obj.anchors <- FindIntegrationAnchors(object.list = list(obj1_1,obj1_2,obj2_1,obj2_2,obj3_1,obj3_2,obj4,obj5,obj6,obj7), dims = 1:20)
obj.combined <- IntegrateData(anchorset = obj.anchors, dims = 1:20)
saveRDS(obj.combined, file="data_integrate.rds")
########
DefaultAssay(obj.combined) <- "integrated"
obj.combined <- FindNeighbors(obj.combined,reduction = "pca",dims = 1:30)
obj.combined <- FindClusters(obj.combined,resolution = 0.6)
obj.combined <- RunTSNE(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:30)
saveRDS(obj.combined, file="data_resolu0.6.rds")
###
pdf("UMAP_0.6.pdf",width=15,height=7)
p1 <- DimPlot(obj.combined, reduction = "umap", group.by="sample",label = FALSE)
p2 <- DimPlot(obj.combined, reduction = "umap",label = TRUE)
CombinePlots(plots  = list(p1, p2))
dev.off()

pdf("UMAP_Sample.pdf",width=15,height=7)
p1 <- DimPlot(obj.combined, reduction = "umap", group.by="group",label = FALSE)
p2 <- DimPlot(obj.combined, reduction = "umap", group.by="Phase",label = FALSE)
CombinePlots(plots  = list(p1, p2))
dev.off()
##marker genes
DefaultAssay(obj.combined) <- "RNA"
obj.combined.markers <- FindAllMarkers(object = obj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(obj.combined.markers,"markers_resolu0.6.txt",row.names=TRUE,col.names=TRUE)
##markers_resolu0.6.txt as Supplementary_File1.xls
