#########
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(infercnv)
## retina as reference
## input retinal dataset
retina <-readRDS("./retinal/data.rds")
current.cluster.ids <- c("0","1","2","4","10","12","3","6","14","9","5","8","7","15","11","13")
new.cluster.ids <- c("Rod","Rod","Rod","Rod","Cone","Cone","Muller","Muller","Astro","microglia","BC","BC","RGC","RGC","AC","HC")
retina@meta.data$type <- plyr::mapvalues(x = retina@meta.data[,"RNA_snn_res.0.6"], from = current.cluster.ids, to = new.cluster.ids)
current.cluster.ids <- c("0","1","2","4","10","12","3","6","14","9","5","8","7","15","11","13")
new.cluster.ids <- c("Rod_C0","Rod_C1","Rod_C2","Rod_c4","Cone_C10","Cone_C12",
"Muller_C3","Muller_C6","Astro","microglia","BC_C5","BC_C8","RGC_C7","RGC_C15","AC","HC")
retina@meta.data$type1 <- plyr::mapvalues(x = retina@meta.data[,"RNA_snn_res.0.6"], from = current.cluster.ids, to = new.cluster.ids)
retina
## input retinoblastoma dataset
obj.combined <-readRDS("./Analysis/data_resolu0.6.rds")
current.cluster.ids <- c("RB1-rep1","RB1-rep2","RB2-rep1","RB2-rep2","RB3-rep1","RB3-rep2","RB4","RB5","RB6","RB7")
new.cluster.ids <- c("RB1","RB1","RB2","RB2","RB3","RB3","RB4","RB5","RB6","RB7")
obj.combined@meta.data$group <- plyr::mapvalues(x = obj.combined@meta.data[,"sample"], from = current.cluster.ids, to = new.cluster.ids)
current.cluster.ids <- c("0","1","2","5","6","3","7","8","9","10","13","16","11","12","15","4","14")
new.cluster.ids <- c("CP-C0","CP-C1","CP-C2","CP-C5","CP-C6","MKI67+ CP-C3","MKI67+ CP-C7","MKI67+ CP-C8","MKI67+ CP-C9","MKI67+ CP-C10","Rod-like-C13","Cone-like-C16","Neural cell-C11","CAF-C12","Glial-C15","Other-C4","Other-C14")
obj.combined@meta.data$celltype <- plyr::mapvalues(x = obj.combined@meta.data[,"integrated_snn_res.0.6"], from = current.cluster.ids, to = new.cluster.ids)
current.cluster.ids <- c("0","1","2","5","6","3","7","8","9","10","13","16","11","12","15","4","14")
new.cluster.ids <- c("CP","CP","CP","CP","CP","MKI67+ CP","MKI67+ CP","MKI67+ CP","MKI67+ CP","MKI67+ CP","Rod-like","Cone-like","Neural cell","CAF","Glial","Other","Other")
obj.combined@meta.data$celltype1 <- plyr::mapvalues(x = obj.combined@meta.data[,"integrated_snn_res.0.6"], from = current.cluster.ids, to = new.cluster.ids)
## merge
pbmc.combined <- merge(obj.combined, y =retina, add.cell.ids = c("RB", "Retina"), project = "Merge")
pbmc.combined@meta.data[69820:69825,]
pbmc.combined@meta.data$Sam <-c(pbmc.combined@meta.data[1:69820,17],pbmc.combined@meta.data[69821:74477,21])
pbmc.combined@meta.data$Newtype <-c(pbmc.combined@meta.data[1:69820,19],pbmc.combined@meta.data[69821:74477,23])
pbmc.combined@meta.data$Newtype1 <-c(pbmc.combined@meta.data[1:69820,18],pbmc.combined@meta.data[69821:74477,24])
setwd("./InferCNV/RB")
saveRDS(pbmc.combined,"data.rds")

##RB
dfcount = as.data.frame(pbmc.combined@assays$RNA@counts) 
groupinfo= data.frame(cellId = colnames(dfcount),cellType= pbmc.combined@meta.data$Newtype1 )# 第二个文件:注释文件，记录肿瘤和正常细胞
library(AnnoProbe)
geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),] # 第三文件基因注释文件
dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
# output
setwd("./InferCNV/RB")
write.table(dfcount ,file ='expFile.txt',sep = '\t',quote = F)
write.table(groupinfo,file = 'metaFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file ='geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)


###run inferCNV
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./InferCNV/RB/expFile.txt",
                                    annotations_file="./InferCNV/RB/metaFiles.txt",
                                    delim="\t",
                                    gene_order_file= "./InferCNV/RB/geneFile.txt",
                                    ref_group_names=c("Retina1","Retina2","Retina3")) # retina as reference 
                                    future::plan("multiprocess",workers=15)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='./InferCNV/RB/output/', 
                             cluster_by_groups=TRUE,
                             denoise=TRUE, 
                             HMM=TRUE, 
							 output_format = "pdf")
		

##Glial as reference
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="./InferCNV/RB/expFile.txt",
                                    annotations_file="./InferCNV/RB/metaFiles.txt",
                                    delim="\t",
                                    gene_order_file= "./InferCNV/RB/geneFile.txt",
                                    ref_group_names=c("Glial-C15")) # Glial as reference
                                    future::plan("multiprocess",workers=15)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='./InferCNV/RB/output/', 
                             cluster_by_groups=TRUE,
                             denoise=TRUE, 
                             HMM=TRUE, 
							 output_format = "pdf")

