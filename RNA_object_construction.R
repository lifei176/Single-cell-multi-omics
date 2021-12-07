setwd("E:/lifei/TPPRC/sc/analysis/Aggregation/")
getwd()
{library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(EnsDb.Mmusculus.v79)
library(Signac)}

#WT 
{inputdata.10x.WT.1 <- Read10X_h5("data/WT_1//filtered_feature_bc_matrix.h5")
inputdata.10x.WT.2 <- Read10X_h5("data/WT_2/filtered_feature_bc_matrix.h5")
rna_counts.WT.1 <- inputdata.10x.WT.1$`Gene Expression`
rna_counts.WT.2 <- inputdata.10x.WT.2$`Gene Expression`
rna.WT.1 <- CreateSeuratObject(counts = rna_counts.WT.1,project = "rna.WT.1")
rna.WT.2 <- CreateSeuratObject(counts = rna_counts.WT.2,project = "rna.WT.2")

#2W 
inputdata.10x.2W.1 <- Read10X_h5("data/2W_1//filtered_feature_bc_matrix.h5")
inputdata.10x.2W.2 <- Read10X_h5("data/2W_2//filtered_feature_bc_matrix.h5")
rna_counts.2W.1 <- inputdata.10x.2W.1$`Gene Expression`
rna_counts.2W.2 <- inputdata.10x.2W.2$`Gene Expression`
rna.2W.1 <- CreateSeuratObject(counts = rna_counts.2W.1,project = "rna.2W.1")
rna.2W.2 <- CreateSeuratObject(counts = rna_counts.2W.2,project = "rna.2W.2")

#1M
inputdata.10x.1M.1 <- Read10X_h5("data/1M_1/filtered_feature_bc_matrix.h5")
inputdata.10x.1M.2 <- Read10X_h5("data/1M_2/filtered_feature_bc_matrix.h5")
rna_counts.1M.1 <- inputdata.10x.1M.1$`Gene Expression`
rna_counts.1M.2 <- inputdata.10x.1M.2$`Gene Expression`
rna.1M.1 <- CreateSeuratObject(counts = rna_counts.1M.1,project = "rna.1M.1")
rna.1M.2 <- CreateSeuratObject(counts = rna_counts.1M.2,project = "rna.1M.2")

#2.5M 
inputdata.10x.2.5M.1 <- Read10X_h5("data/2_5M_1/filtered_feature_bc_matrix.h5")
inputdata.10x.2.5M.2 <- Read10X_h5("data/2_5M_2/filtered_feature_bc_matrix.h5")
inputdata.10x.2.5M.3 <- Read10X_h5("data/2_5M_3/filtered_feature_bc_matrix.h5")
inputdata.10x.2.5M.4 <- Read10X_h5("data/2_5M_4/filtered_feature_bc_matrix.h5")
rna_counts.2.5M.1 <- inputdata.10x.2.5M.1$`Gene Expression`
rna_counts.2.5M.2 <- inputdata.10x.2.5M.2$`Gene Expression`
rna_counts.2.5M.3 <- inputdata.10x.2.5M.3$`Gene Expression`
rna_counts.2.5M.4 <- inputdata.10x.2.5M.4$`Gene Expression`
rna.2.5M.1 <- CreateSeuratObject(counts = rna_counts.2.5M.1,project = "rna.2.5M.1")
rna.2.5M.2 <- CreateSeuratObject(counts = rna_counts.2.5M.2,project = "rna.2.5M.2")
rna.2.5M.3 <- CreateSeuratObject(counts = rna_counts.2.5M.3,project = "rna.2.5M.3")
rna.2.5M.4 <- CreateSeuratObject(counts = rna_counts.2.5M.4,project = "rna.2.5M.4")

#3.5M 
inputdata.10x.3.5M.1 <- Read10X_h5("data/3_5M_1/filtered_feature_bc_matrix.h5")
inputdata.10x.3.5M.2 <- Read10X_h5("data/3_5M_2/filtered_feature_bc_matrix.h5")
rna_counts.3.5M.1 <- inputdata.10x.3.5M.1$`Gene Expression`
rna_counts.3.5M.2 <- inputdata.10x.3.5M.2$`Gene Expression`
rna.3.5M.1 <- CreateSeuratObject(counts = rna_counts.3.5M.1,project = "rna.3.5M.1")
rna.3.5M.2 <- CreateSeuratObject(counts = rna_counts.3.5M.2,project = "rna.3.5M.2")

#4.5M 
inputdata.10x.4.5M.1 <- Read10X_h5("data/4_5M_1/filtered_feature_bc_matrix.h5")
inputdata.10x.4.5M.2 <- Read10X_h5("data/4_5M_2/filtered_feature_bc_matrix.h5")
rna_counts.4.5M.1 <- inputdata.10x.4.5M.1$`Gene Expression`
rna_counts.4.5M.2 <- inputdata.10x.4.5M.2$`Gene Expression`
rna.4.5M.1 <- CreateSeuratObject(counts = rna_counts.4.5M.1,project = "rna.4.5M.1")
rna.4.5M.2 <- CreateSeuratObject(counts = rna_counts.4.5M.2,project = "rna.4.5M.2")

#RC
inputdata.10x.RC.1 <- Read10X_h5("data/RC/filtered_feature_bc_matrix.h5")
rna_counts.RC.1 <- inputdata.10x.RC.1$`Gene Expression`
rna.RC <- CreateSeuratObject(counts = rna_counts.RC.1,project = "rna.RC")
# Create Seurat object

###Merge Seurat object and remove data to release storage
rna.All.combined <- merge(rna.WT.1, y = c(rna.WT.2,rna.2W.1,rna.2W.2,rna.1M.1,rna.1M.2,rna.2.5M.1,rna.2.5M.2,rna.2.5M.3,rna.2.5M.4,rna.3.5M.1,rna.3.5M.2,rna.4.5M.1,rna.4.5M.2,rna.RC), add.cell.ids = c("rna.WT.1", "rna.WT.2","rna.2W.1","rna.2W.2","rna.1M.1","rna.1M.2","rna.2.5M.1","rna.2.5M.2","rna.2.5M.3","rna.2.5M.4","rna.3.5M.1","rna.3.5M.2","rna.4.5M.1","rna.4.5M.2","rna.RC"), project = "rna.All.combined")
}
rm(inputdata.10x.WT.1,inputdata.10x.WT.2,inputdata.10x.2W.1,inputdata.10x.2W.2)
rm(inputdata.10x.1M.1,inputdata.10x.1M.2,inputdata.10x.2.5M.1,inputdata.10x.2.5M.2,inputdata.10x.2.5M.3,inputdata.10x.2.5M.4)
rm(inputdata.10x.3.5M.1,inputdata.10x.3.5M.2,inputdata.10x.4.5M.1,inputdata.10x.4.5M.2,inputdata.10x.RC.1)
rm(rna.WT.1,rna.WT.2,rna.2W.1,rna.2W.2,rna.1M.1,rna.1M.2,rna.2.5M.1,rna.2.5M.2,rna.2.5M.3,rna.2.5M.4,rna.3.5M.1,rna.3.5M.2,rna.4.5M.1,rna.4.5M.2,rna.RC)
rm(rna_counts.WT.1,rna_counts.WT.2,rna_counts.2W.1,rna_counts.2W.2,rna_counts.1M.1,rna_counts.1M.2,rna_counts.2.5M.1,rna_counts.2.5M.2)
rm(rna_counts.2.5M.3,rna_counts.2.5M.4,rna_counts.3.5M.1,rna_counts.3.5M.2,rna_counts.4.5M.1,rna_counts.4.5M.2,rna_counts.RC.1)
gc()
#unique(sapply(X = strsplit(colnames(rna.All.combined), split = "_"), FUN = "[", 1))
##mt gene marking
table(rna.All.combined@meta.data$orig.ident)
rna.All.combined[["percent.mt"]] <- PercentageFeatureSet(rna.All.combined, pattern = "^mt-")

### QC
VlnPlot(rna.All.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
{plot1 <- FeatureScatter(rna.All.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(rna.All.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2}
###
rna.All.combined.list <- SplitObject(rna.All.combined, split.by = "orig.ident")
  rna.All.combined.list <- lapply(X = rna.All.combined.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
  })
  
  features <- SelectIntegrationFeatures(object.list = rna.All.combined.list)
  rna.All.combined.list <- lapply(X = rna.All.combined.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
#######remove rna.All.combined 
  rm(rna.All.combined,anchors)
  gc()
  # reference = c(1, 2)would reduce time
  #anchors <- FindIntegrationAnchors(object.list = rna.All.combined.list, reference = c(1, 2), reduction = "rpca",dims = 1:50) 
                                  
  #all pairwise anchors
  anchors <- FindIntegrationAnchors(object.list = rna.All.combined.list, reduction = "rpca", 
                                    dims = 1:20)

  #This step would waste much time
  rna.All.combined.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
    rm(rna.All.combined.list,)
  gc()
  rna.All.combined.integrated <- ScaleData(rna.All.combined.integrated, verbose = FALSE)
  rna.All.combined.integrated <- RunPCA(rna.All.combined.integrated, verbose = FALSE)
  rna.All.combined.integrated <- RunUMAP(rna.All.combined.integrated, dims = 1:50)
  rna.All.combined.integrated <- FindNeighbors(rna.All.combined.integrated, dims = 1:50)
  rna.All.combined.integrated <- FindClusters(rna.All.combined.integrated, resolution = 0.1)#15 clusters

#DimPlot(rna.All.combined.integrated, group.by = "orig.ident",split.by="orig.ident")
#DimPlot(rna.All.combined.integrated, group.by = "orig.ident")
write.csv(rna.All.combined.integrated@meta.data,"rna.All.combined.integrated_meta.data.csv")  
###definition clusters
timepoint<-c(rep("1_WT",times=18288),rep("2_2W",times=8479),rep("3_1M",times=15800),rep("4_2.5M",times=48841),rep("5_3.5M",times=25346),rep("6_4.5M",times=31104),rep("7_RC",times=8220))    
rna.All.combined.integrated$timepoint <-timepoint  

#####Low_quality_cell identification
rna.All.combined.integrated <- RenameIdents(rna.All.combined.integrated, '0' = 'Low_quality_cell','1' = 'Luminal','2' = 'Neuroendocrine','3' = 'Stromal_1','4' = 'Neutrophil')
rna.All.combined.integrated <- RenameIdents(rna.All.combined.integrated ,'5' = 'Macrophage', '6' = 'Basal','7' = 'Endothelial cell','8' = 'T cell','9' = 'Seminal vesicle','10' = 'Adipocyte','11' = 'Stromal_2','12' = 'B cell','13' = 'Reln_high','14' = 'Neuron')

#####################################
#####################################
#####################################subset to remove low-quality cells
rna.All.combined.integrated.filter<-subset(rna.All.combined.integrated,integrated_snn_res.0.1!="0")
rna.All.combined.integrated.filter#107201
head(rna.All.combined.integrated.filter@meta.data)
rm(rna.All.combined.integrated)
{rna.All.combined.integrated.filter <- ScaleData(rna.All.combined.integrated.filter, verbose = FALSE)
rna.All.combined.integrated.filter <- RunPCA(rna.All.combined.integrated.filter, verbose = FALSE)
rna.All.combined.integrated.filter <- RunUMAP(rna.All.combined.integrated.filter, dims = 1:50)
rna.All.combined.integrated.filter <- FindNeighbors(rna.All.combined.integrated.filter, dims = 1:50)
rna.All.combined.integrated.filter <- FindClusters(rna.All.combined.integrated.filter, resolution = 0.055)}#13 clusters
save(rna.All.combined.integrated.filter,file="rna.All.combined.integrated.filter.Rdata")



























