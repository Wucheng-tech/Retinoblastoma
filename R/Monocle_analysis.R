#########monocle new choose
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(HSMMSingleCell)
library(monocle)
library(Seurat)
## input retinoblastoma retinal
obj.combined <-readRDS("./retinal/intergrate_RB/data_integrate1.rds")
ind1 <-which(obj.combined@meta.data[,17]=="CP")
ind <-sample(ind1,length(ind1)*0.1)
ind1 <-which(obj.combined@meta.data[,17]=="MKI67+ CP")
index <-sample(ind1,length(ind1)*0.05)
index2<-NULL
for(i in c("Cone","Cone-like")){
ind1 <-which(obj.combined@meta.data[,17]==i)
ind2 <-sample(ind1,length(ind1))
index2 <-c(ind2,index2)
}
index3 <-c(ind,index,index2)

pbmc <-obj.combined
data <- as(as.matrix(pbmc@assays$RNA@counts[,index3]), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = pbmc@meta.data[index3,])
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
HSMM<-monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 500))
print(head(pData(HSMM)))

###choosing genes that define progress
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~Newtype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)
setwd("./monocle")
HSMM <-readRDS("RB_Newtype_sub.rds")
p1 <- plot_cell_trajectory(HSMM, color_by = "Newtype")
p2 <- plot_cell_trajectory(HSMM, color_by = "Newtype1")
pdf("pseudotime_RB_Newtype_sub.pdf",width=15,height=7)
CombinePlots(plots = list(p1, p2))
dev.off()
saveRDS(HSMM, file="RB_Newtype_sub.rds")
pdf("pseudotime_RB_wrap_sub.pdf",width=15,height=8)
plot_cell_trajectory(HSMM, color_by = "Newtype") +facet_wrap(~Newtype, nrow = 3)
dev.off()
pdf("pseudotime_RB_wrap_sub1.pdf",width=15,height=8)
plot_cell_trajectory(HSMM, color_by = "Newtype1") +facet_wrap(~Newtype1, nrow = 3)
dev.off()
pdf("pseudotime_RB_state_sub.pdf",width=15,height=8)
plot_cell_trajectory(HSMM, color_by = "State") +facet_wrap(~State, nrow = 3)
dev.off()

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Newtype1)[,c(1)]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM <- orderCells(HSMM, root_state = GM_state(HSMM))
pdf("pseudotime_GM_sub.pdf",width=15,height=8)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()
pdf("pseudotime_GM_state.pdf",width=15,height=8)
plot_cell_trajectory(HSMM, color_by = "State") +facet_wrap(~State, nrow = 1)
dev.off()
p1 <- plot_cell_trajectory(HSMM, color_by = "Newtype1")
p2 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
pdf("pseudotime_GM_sub1.pdf",width=18,height=7)
CombinePlots(plots = list(p1, p2))
dev.off()

blast_genes <- row.names(subset(fData(HSMM),gene_short_name %in% c("MKI67","KIF14","TOP2A","CDC20","CDCA8")))
pdf("pseudotime_Gene.pdf",width=10,height=5)
plot_genes_jitter(HSMM[blast_genes,],grouping = "State", min_expr = 0.1)
dev.off()

HSMM_expressed_genes <-  row.names(subset(fData(HSMM),num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),gene_short_name %in% c("MKI67","TOP2A")))
cds_subset <- HSMM_filtered[my_genes,]
pdf("pseudotime_Gene1.pdf",width=9,height=5)
plot_genes_in_pseudotime(cds_subset, color_by = "Newtype1")
dev.off()

HSMM_expressed_genes <-  row.names(subset(fData(HSMM),num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),gene_short_name %in% c("PDE6H","GNGT2")))
cds_subset <- HSMM_filtered[my_genes,]
pdf("pseudotime_Gene2.pdf",width=9,height=5)
plot_genes_in_pseudotime(cds_subset, color_by = "Newtype1")
dev.off()
