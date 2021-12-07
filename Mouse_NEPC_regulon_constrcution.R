library(SCENIC)
library(foreach)
library(Seurat)
library(arrow)
###############Sampling 0.2
###############Sampling 0.2
###############Sampling 0.2
###############Sampling 0.2
###############Sampling 0.2
###############Sampling 0.2
setwd("E:\lifei\TPPRC\sc\analysis\Aggregation\outs\Linux\Integrated\20211020\Mouse_regulon_reanalysis")
load("/sibcb2/gaodonglab2/lifei/NEPC/sc/R_analysis/rna.All.combined.integrated.filter.prostate.filter.20210604.Rdata")
######Downsamping 500 cells for each cluster
seu<-subset(rna.All.combined.integrated.filter.prostate.filter,downsample=500)
exprMat<-seu@assays$RNA@data[which(rowSums(seu@assays$RNA@counts) > 1000),]
dim(exprMat)
exprMat<-as.matrix(exprMat)
### cell information
cellinfo<-data.frame(rna.All.combined.integrated.filter.prostate.filter@meta.data)
cellinfo<-cellinfo[,c("orig.ident","seurat_clusters")]
colnames(cellinfo)<-c("sample","cell_type")
saveRDS(cellinfo,file="int/cellinfo.Rds")
#
scenicOptions <- initializeScenic(org="hgnc", dbDir="/sibcb2/gaodonglab2/lifei/NEPC/sc/R_analysis/Human_NEPC_regulon/cisTarget_databases", nCores=6)
scenicOptions@inputDatasetInfo$cellinfo<-"int/cellinfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions,minCountsPerGene = 3*0.01*ncol(exprMat),
                           minSamples = ncol(exprMat)*0.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
####correlation calculation
runCorrelation(exprMat_filtered, scenicOptions)
####TF-targets regression
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
### Build and score the GRN
#scenicOptions=readRDS(file="int/scenicOptions.Rds")
#scenicOptions@settings$verbose <- TRUE
#scenicOptions@settings$nCores <-1
#scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

