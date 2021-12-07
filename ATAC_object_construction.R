setwd("E:/lifei/TPPRC/sc/analysis/Aggregation/outs/Figure/Integrative_ATAC/")
getwd()
{library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(EnsDb.Mmusculus.v79)
  library(Signac)
  library(S4Vectors)
  library(patchwork)
  set.seed(1234)}

#WT 
{inputdata.10x.WT.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/WT_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.WT.2 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/WT_2/filtered_feature_bc_matrix.h5")
  atac_counts.WT.1 <- inputdata.10x.WT.1$Peaks
  atac_counts.WT.2 <- inputdata.10x.WT.2$Peaks

  #2W 
  inputdata.10x.2W.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/2W_1//filtered_feature_bc_matrix.h5")
  inputdata.10x.2W.2 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/2W_2//filtered_feature_bc_matrix.h5")
  atac_counts.2W.1 <- inputdata.10x.2W.1$Peaks
  atac_counts.2W.2 <- inputdata.10x.2W.2$Peaks

  #1M
  inputdata.10x.1M.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/1M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.1M.2 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/1M_2/filtered_feature_bc_matrix.h5")
  atac_counts.1M.1 <- inputdata.10x.1M.1$Peaks
  atac_counts.1M.2 <- inputdata.10x.1M.2$Peaks

  #2.5M 
  inputdata.10x.2.5M.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.2.5M.2 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_2/filtered_feature_bc_matrix.h5")
  inputdata.10x.2.5M.3 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_3/filtered_feature_bc_matrix.h5")
  inputdata.10x.2.5M.4 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_4/filtered_feature_bc_matrix.h5")
  atac_counts.2.5M.1 <- inputdata.10x.2.5M.1$Peaks
  atac_counts.2.5M.2 <- inputdata.10x.2.5M.2$Peaks
  atac_counts.2.5M.3 <- inputdata.10x.2.5M.3$Peaks
  atac_counts.2.5M.4 <- inputdata.10x.2.5M.4$Peaks

  #3.5M 
  inputdata.10x.3.5M.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/3_5M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.3.5M.2 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/3_5M_2/filtered_feature_bc_matrix.h5")
  atac_counts.3.5M.1 <- inputdata.10x.3.5M.1$Peaks
  atac_counts.3.5M.2 <- inputdata.10x.3.5M.2$Peaks

  #4.5M 
  inputdata.10x.4.5M.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/4_5M_1/filtered_feature_bc_matrix.h5")
  inputdata.10x.4.5M.2 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/4_5M_2/filtered_feature_bc_matrix.h5")
  atac_counts.4.5M.1 <- inputdata.10x.4.5M.1$Peaks
  atac_counts.4.5M.2 <- inputdata.10x.4.5M.2$Peaks

  #RC
  inputdata.10x.RC.1 <- Read10X_h5("E:/lifei/TPPRC/sc/analysis/Aggregation/data/RC/filtered_feature_bc_matrix.h5")
  atac_counts.RC.1 <- inputdata.10x.RC.1$Peaks}

#######Define frag.files
{frag.file.WT.1 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/WT_1/outs/atac_fragments.tsv.gz"
  frag.file.WT.2 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/WT_2/outs/atac_fragments.tsv.gz"
  frag.file.2W.1 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/2W_1/outs/atac_fragments.tsv.gz"
  frag.file.2W.2 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/2W_2/outs/atac_fragments.tsv.gz"
  frag.file.1M.1 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/1M_1/outs/atac_fragments.tsv.gz"
  frag.file.1M.2 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/1M_2/outs/atac_fragments.tsv.gz"
  frag.file.2.5M.1 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/2_5M_1/outs/atac_fragments.tsv.gz"
  frag.file.2.5M.2 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/2_5M_2/outs/atac_fragments.tsv.gz"
  frag.file.2.5M.3 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/2_5M_3/outs/atac_fragments.tsv.gz"
  frag.file.2.5M.4 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/2_5M_4/outs/atac_fragments.tsv.gz"
  frag.file.3.5M.1 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/3_5M_1/outs/atac_fragments.tsv.gz"
  frag.file.3.5M.2 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/3_5M_2/outs/atac_fragments.tsv.gz"
  frag.file.4.5M.1 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/4_5M_1/outs/atac_fragments.tsv.gz"
  frag.file.4.5M.2 <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/4_5M_2/outs/atac_fragments.tsv.gz"
  frag.file.RC <- "/sibcb2/gaodonglab2/lifei/NEPC/raw_data/SC_RNA/log/RC/outs/atac_fragments.tsv.gz"}
#Create Seurat object
#construct atac seurat object
#atac.WT.1
{{grange.counts.WT.1 <- StringToGRanges(rownames(atac_counts.WT.1), sep = c(":", "-"))
grange.use.WT.1 <- seqnames(grange.counts.WT.1) %in% standardChromosomes(grange.counts.WT.1)
atac_counts.WT.1 <- atac_counts.WT.1[as.vector(grange.use.WT.1), ]
annotations.WT.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations.WT.1) <- 'UCSC'
genome(annotations.WT.1) <- "mm10"

chrom_assay.WT.1 <- CreateChromatinAssay(
  counts = atac_counts.WT.1,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file.WT.1,
  min.cells = 10,
  annotation = annotations.WT.1
)}
  #atac.WT.2
  {grange.counts.WT.2 <- StringToGRanges(rownames(atac_counts.WT.2), sep = c(":", "-"))
    grange.use.WT.2 <- seqnames(grange.counts.WT.2) %in% standardChromosomes(grange.counts.WT.2)
    atac_counts.WT.2 <- atac_counts.WT.2[as.vector(grange.use.WT.2), ]
    annotations.WT.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.WT.2) <- 'UCSC'
    genome(annotations.WT.2) <- "mm10"
    
    chrom_assay.WT.2 <- CreateChromatinAssay(
      counts = atac_counts.WT.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.WT.2,
      min.cells = 10,
      annotation = annotations.WT.2
    )}
  
  #atac.2W.1
  {grange.counts.2W.1 <- StringToGRanges(rownames(atac_counts.2W.1), sep = c(":", "-"))
    grange.use.2W.1 <- seqnames(grange.counts.2W.1) %in% standardChromosomes(grange.counts.2W.1)
    atac_counts.2W.1 <- atac_counts.2W.1[as.vector(grange.use.2W.1), ]
    annotations.2W.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2W.1) <- 'UCSC'
    genome(annotations.2W.1) <- "mm10"
    
    chrom_assay.2W.1 <- CreateChromatinAssay(
      counts = atac_counts.2W.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2W.1,
      min.cells = 10,
      annotation = annotations.2W.1
    )}
  
  
  #atac.2W.2
  {grange.counts.2W.2 <- StringToGRanges(rownames(atac_counts.2W.2), sep = c(":", "-"))
    grange.use.2W.2 <- seqnames(grange.counts.2W.2) %in% standardChromosomes(grange.counts.2W.2)
    atac_counts.2W.2 <- atac_counts.2W.2[as.vector(grange.use.2W.2), ]
    annotations.2W.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2W.2) <- 'UCSC'
    genome(annotations.2W.2) <- "mm10"
    
    chrom_assay.2W.2 <- CreateChromatinAssay(
      counts = atac_counts.2W.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2W.2,
      min.cells = 10,
      annotation = annotations.2W.2
    )}
  
  #atac.1M.1
  {grange.counts.1M.1 <- StringToGRanges(rownames(atac_counts.1M.1), sep = c(":", "-"))
    grange.use.1M.1 <- seqnames(grange.counts.1M.1) %in% standardChromosomes(grange.counts.1M.1)
    atac_counts.1M.1 <- atac_counts.1M.1[as.vector(grange.use.1M.1), ]
    annotations.1M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.1M.1) <- 'UCSC'
    genome(annotations.1M.1) <- "mm10"
    
    chrom_assay.1M.1 <- CreateChromatinAssay(
      counts = atac_counts.1M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.1M.1,
      min.cells = 10,
      annotation = annotations.1M.1
    )}
  
  #atac.1M.2
  {grange.counts.1M.2 <- StringToGRanges(rownames(atac_counts.1M.2), sep = c(":", "-"))
    grange.use.1M.2 <- seqnames(grange.counts.1M.2) %in% standardChromosomes(grange.counts.1M.2)
    atac_counts.1M.2 <- atac_counts.1M.2[as.vector(grange.use.1M.2), ]
    annotations.1M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.1M.2) <- 'UCSC'
    genome(annotations.1M.2) <- "mm10"
    
    chrom_assay.1M.2 <- CreateChromatinAssay(
      counts = atac_counts.1M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.1M.2,
      min.cells = 10,
      annotation = annotations.1M.2
    )}
  #atac.2.5M.1
  {grange.counts.2.5M.1 <- StringToGRanges(rownames(atac_counts.2.5M.1), sep = c(":", "-"))
    grange.use.2.5M.1 <- seqnames(grange.counts.2.5M.1) %in% standardChromosomes(grange.counts.2.5M.1)
    atac_counts.2.5M.1 <- atac_counts.2.5M.1[as.vector(grange.use.2.5M.1), ]
    annotations.2.5M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.1) <- 'UCSC'
    genome(annotations.2.5M.1) <- "mm10"
    
    chrom_assay.2.5M.1 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.1,
      min.cells = 10,
      annotation = annotations.2.5M.1
    )}
  
  #atac.2.5M.2
  {grange.counts.2.5M.2 <- StringToGRanges(rownames(atac_counts.2.5M.2), sep = c(":", "-"))
    grange.use.2.5M.2 <- seqnames(grange.counts.2.5M.2) %in% standardChromosomes(grange.counts.2.5M.2)
    atac_counts.2.5M.2 <- atac_counts.2.5M.2[as.vector(grange.use.2.5M.2), ]
    annotations.2.5M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.2) <- 'UCSC'
    genome(annotations.2.5M.2) <- "mm10"
    
    chrom_assay.2.5M.2 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.2,
      min.cells = 10,
      annotation = annotations.2.5M.2
    )}
  
  #atac.2.5M.3
  {grange.counts.2.5M.3 <- StringToGRanges(rownames(atac_counts.2.5M.3), sep = c(":", "-"))
    grange.use.2.5M.3 <- seqnames(grange.counts.2.5M.3) %in% standardChromosomes(grange.counts.2.5M.3)
    atac_counts.2.5M.3 <- atac_counts.2.5M.3[as.vector(grange.use.2.5M.3), ]
    annotations.2.5M.3 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.3) <- 'UCSC'
    genome(annotations.2.5M.3) <- "mm10"
    
    chrom_assay.2.5M.3 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.3,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.3,
      min.cells = 10,
      annotation = annotations.2.5M.3
    )}
  #atac.2.5M.4
  {grange.counts.2.5M.4 <- StringToGRanges(rownames(atac_counts.2.5M.4), sep = c(":", "-"))
    grange.use.2.5M.4 <- seqnames(grange.counts.2.5M.4) %in% standardChromosomes(grange.counts.2.5M.4)
    atac_counts.2.5M.4 <- atac_counts.2.5M.4[as.vector(grange.use.2.5M.4), ]
    annotations.2.5M.4 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.2.5M.4) <- 'UCSC'
    genome(annotations.2.5M.4) <- "mm10"
    
    chrom_assay.2.5M.4 <- CreateChromatinAssay(
      counts = atac_counts.2.5M.4,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.2.5M.4,
      min.cells = 10,
      annotation = annotations.2.5M.4
    )}
  
  #atac.3.5M.1
  {grange.counts.3.5M.1 <- StringToGRanges(rownames(atac_counts.3.5M.1), sep = c(":", "-"))
    grange.use.3.5M.1 <- seqnames(grange.counts.3.5M.1) %in% standardChromosomes(grange.counts.3.5M.1)
    atac_counts.3.5M.1 <- atac_counts.3.5M.1[as.vector(grange.use.3.5M.1), ]
    annotations.3.5M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.3.5M.1) <- 'UCSC'
    genome(annotations.3.5M.1) <- "mm10"
    
    chrom_assay.3.5M.1 <- CreateChromatinAssay(
      counts = atac_counts.3.5M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.3.5M.1,
      min.cells = 10,
      annotation = annotations.3.5M.1
    )}
  
  #atac.3.5M.2
  {grange.counts.3.5M.2 <- StringToGRanges(rownames(atac_counts.3.5M.2), sep = c(":", "-"))
    grange.use.3.5M.2 <- seqnames(grange.counts.3.5M.2) %in% standardChromosomes(grange.counts.3.5M.2)
    atac_counts.3.5M.2 <- atac_counts.3.5M.2[as.vector(grange.use.3.5M.2), ]
    annotations.3.5M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.3.5M.2) <- 'UCSC'
    genome(annotations.3.5M.2) <- "mm10"
    
    chrom_assay.3.5M.2 <- CreateChromatinAssay(
      counts = atac_counts.3.5M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.3.5M.2,
      min.cells = 10,
      annotation = annotations.3.5M.2
    )}
  
  #atac.4.5M.1
  {grange.counts.4.5M.1 <- StringToGRanges(rownames(atac_counts.4.5M.1), sep = c(":", "-"))
    grange.use.4.5M.1 <- seqnames(grange.counts.4.5M.1) %in% standardChromosomes(grange.counts.4.5M.1)
    atac_counts.4.5M.1 <- atac_counts.4.5M.1[as.vector(grange.use.4.5M.1), ]
    annotations.4.5M.1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.4.5M.1) <- 'UCSC'
    genome(annotations.4.5M.1) <- "mm10"
    
    chrom_assay.4.5M.1 <- CreateChromatinAssay(
      counts = atac_counts.4.5M.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.4.5M.1,
      min.cells = 10,
      annotation = annotations.4.5M.1
    )}
  
  #atac.4.5M.2
  {grange.counts.4.5M.2 <- StringToGRanges(rownames(atac_counts.4.5M.2), sep = c(":", "-"))
    grange.use.4.5M.2 <- seqnames(grange.counts.4.5M.2) %in% standardChromosomes(grange.counts.4.5M.2)
    atac_counts.4.5M.2 <- atac_counts.4.5M.2[as.vector(grange.use.4.5M.2), ]
    annotations.4.5M.2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.4.5M.2) <- 'UCSC'
    genome(annotations.4.5M.2) <- "mm10"
    
    chrom_assay.4.5M.2 <- CreateChromatinAssay(
      counts = atac_counts.4.5M.2,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.4.5M.2,
      min.cells = 10,
      annotation = annotations.4.5M.2
    )}
  
  #atac.RC
  {grange.counts.RC <- StringToGRanges(rownames(atac_counts.RC.1), sep = c(":", "-"))
    grange.use.RC <- seqnames(grange.counts.RC) %in% standardChromosomes(grange.counts.RC)
    atac_counts.RC.1 <- atac_counts.RC.1[as.vector(grange.use.RC), ]
    annotations.RC <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations.RC) <- 'UCSC'
    genome(annotations.RC) <- "mm10"
    
    chrom_assay.RC <- CreateChromatinAssay(
      counts = atac_counts.RC.1,
      sep = c(":", "-"),
      genome = 'mm10',
      fragments = frag.file.RC,
      min.cells = 10,
      annotation = annotations.RC
    )}}
####release storage 
rm(atac_counts.1M.1,atac_counts.1M.2,atac_counts.2.5M.1,atac_counts.2.5M.2,atac_counts.2.5M.3,atac_counts.2.5M.4)
rm(atac_counts.2W.1,atac_counts.2W.2,atac_counts.3.5M.1,atac_counts.3.5M.2,atac_counts.4.5M.1,atac_counts.4.5M.2)
rm(atac_counts.WT.1,atac_counts.WT.2,atac_counts.RC.1)
gc()

##############construct consensus_peaksets
#BiocManager::install("DiffBind")
#library(amap)
library(DiffBind)
data<-"SampleSheet.csv"
dbObj <- dba(sampleSheet=data)
consensus_peakset <- dba.peakset(dbObj,dbObj$masks$Consensus,bRetrieve=TRUE)
dbObj_consensus1 <- dba.peakset(dbObj, consensus = DBA_CONDITION, minOverlap=0.5)#
dim(dbObj_consensus1)
consensus_peakset <- dba.peakset(dbObj_consensus1,dbObj_consensus1$masks$Consensus,bRetrieve=TRUE)
consensus_peakset<-consensus_peakset[1:105020,]
length(consensus_peakset)
class(consensus_peakset)
write.csv(consensus_peakset,"peakset_minOverlap_0.5_remove.csv")
##############re-construct peak-counts matrix
{peaks_WT.2 <- FeatureMatrix(
  fragments = Fragments(chrom_assay.WT.2),
  features = consensus_peakset,
  cells = colnames(chrom_assay.WT.2)
)
  
  grange.counts.peakset.inter <- StringToGRanges(rownames(peaks_WT.2), sep = c(":", "-"))
  grange.use.peakset.inter <- seqnames(grange.counts.peakset.inter) %in% standardChromosomes(grange.counts.peakset.inter)
  atac_counts.peakset.inter <- peaks_WT.2[as.vector(grange.use.peakset.inter), ]
  annotations.peakset.inter <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations.peakset.inter) <- 'UCSC'
  genome(annotations.peakset.inter) <- "mm10"
  
  integarted.atac.WT.2.counts<- CreateChromatinAssay(
    counts = peaks_WT.2,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.WT.2
  )
  
  peaks_WT.1 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.WT.1),
    features = consensus_peakset,
    cells = colnames(chrom_assay.WT.1)
  )
  
  integarted.atac.WT.1.counts<- CreateChromatinAssay(
    counts = peaks_WT.1,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.WT.1
  )
  
  peaks_2W.1 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.2W.1),
    features = consensus_peakset,
    cells = colnames(chrom_assay.2W.1)
  )
  integarted.atac.2W.1.counts<- CreateChromatinAssay(
    counts = peaks_2W.1,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.2W.1
  )
  
  
  peaks_2W.2 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.2W.2),
    features = consensus_peakset,
    cells = colnames(chrom_assay.2W.2)
  )
  integarted.atac.2W.2.counts<- CreateChromatinAssay(
    counts = peaks_2W.2,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.2W.2
  )
  
  peaks_1M.1 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.1M.1),
    features = consensus_peakset,
    cells = colnames(chrom_assay.1M.1)
  )
  integarted.atac.1M.1.counts<- CreateChromatinAssay(
    counts = peaks_1M.1,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.1M.1
  )
  
  peaks_1M.2 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.1M.2),
    features = consensus_peakset,
    cells = colnames(chrom_assay.1M.2)
  )
  integarted.atac.1M.2.counts<- CreateChromatinAssay(
    counts = peaks_1M.2,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.1M.2
  )
  
  peaks_2.5M.1 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.2.5M.1),
    features = consensus_peakset,
    cells = colnames(chrom_assay.2.5M.1)
  )
  integarted.atac.2.5M.1.counts<- CreateChromatinAssay(
    counts = peaks_2.5M.1,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.2.5M.1
  )
  peaks_2.5M.2 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.2.5M.2),
    features = consensus_peakset,
    cells = colnames(chrom_assay.2.5M.2)
  )
  integarted.atac.2.5M.2.counts<- CreateChromatinAssay(
    counts = peaks_2.5M.2,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.2.5M.2
  )
  
  peaks_2.5M.3 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.2.5M.3),
    features = consensus_peakset,
    cells = colnames(chrom_assay.2.5M.3)
  )
  integarted.atac.2.5M.3.counts<- CreateChromatinAssay(
    counts = peaks_2.5M.3,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.2.5M.3
  )
  
  peaks_2.5M.4 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.2.5M.4),
    features = consensus_peakset,
    cells = colnames(chrom_assay.2.5M.4)
  )
  integarted.atac.2.5M.4.counts<- CreateChromatinAssay(
    counts = peaks_2.5M.4,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.2.5M.4
  )
  
  peaks_3.5M.1 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.3.5M.1),
    features = consensus_peakset,
    cells = colnames(chrom_assay.3.5M.1)
  )
  integarted.atac.3.5M.1.counts<- CreateChromatinAssay(
    counts = peaks_3.5M.1,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.3.5M.1
  )
  
  peaks_3.5M.2 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.3.5M.2),
    features = consensus_peakset,
    cells = colnames(chrom_assay.3.5M.2)
  )
  integarted.atac.3.5M.2.counts<- CreateChromatinAssay(
    counts = peaks_3.5M.2,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.3.5M.2
  )
  peaks_4.5M.1 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.4.5M.1),
    features = consensus_peakset,
    cells = colnames(chrom_assay.4.5M.1)
  )
  integarted.atac.4.5M.1.counts<- CreateChromatinAssay(
    counts = peaks_4.5M.1,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.4.5M.1
  )
  peaks_4.5M.2 <- FeatureMatrix(
    fragments = Fragments(chrom_assay.4.5M.2),
    features = consensus_peakset,
    cells = colnames(chrom_assay.4.5M.2)
  )
  integarted.atac.4.5M.2.counts<- CreateChromatinAssay(
    counts = peaks_4.5M.2,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.4.5M.2
  )
  peaks_RC <- FeatureMatrix(
    fragments = Fragments(chrom_assay.RC),
    features = consensus_peakset,
    cells = colnames(chrom_assay.RC)
  )
  
  integarted.atac.RC.counts<- CreateChromatinAssay(
    counts = peaks_RC,
    sep = c(":", "-"),
    genome = 'mm10',
    min.cells = 1,
    annotation=annotations.peakset.inter,
    fragments = frag.file.RC
  )}
#memory release 
rm(peaks_1M.1,peaks_1M.2,peaks_2.5M.1,peaks_2.5M.2,peaks_2.5M.3,peaks_2.5M.4,peaks_2W.1,peaks_2W.2,peaks_3.5M.1,peaks_3.5M.2,peaks_4.5M.1,peaks_4.5M.2,peaks_RC,peaks_WT.1,peaks_WT.2)
gc()
rm(inputdata.10x.1M.1,inputdata.10x.1M.2,inputdata.10x.2.5M.1,inputdata.10x.2.5M.2,inputdata.10x.2.5M.3,inputdata.10x.2.5M.4,inputdata.10x.2W.1,inputdata.10x.2W.2,inputdata.10x.3.5M.1,inputdata.10x.3.5M.2,inputdata.10x.4.5M.1,inputdata.10x.4.5M.2,inputdata.10x.RC.1,inputdata.10x.WT.1,inputdata.10x.WT.2)
gc()
rm(chrom_assay.1M.1,chrom_assay.1M.2,chrom_assay.2.5M.1,chrom_assay.2.5M.2,chrom_assay.2.5M.3,chrom_assay.2.5M.4,chrom_assay.2W.1,chrom_assay.2W.2,chrom_assay.3.5M.1,chrom_assay.3.5M.2,chrom_assay.4.5M.1,chrom_assay.4.5M.2,chrom_assay.RC,chrom_assay.WT.1,chrom_assay.WT.2)
gc()
rm(annotations.1M.1,annotations.1M.2,annotations.2.5M.1,annotations.2.5M.2,annotations.2.5M.3,annotations.2.5M.4,annotations.2W.1,annotations.2W.2,annotations.3.5M.1,annotations.3.5M.2,annotations.4.5M.1,annotations.4.5M.2,annotations.RC,annotations.WT.1,annotations.WT.2)
gc()

##########
####1M.1 renames
{integarted.atac.1M.1.renames<-merge(integarted.atac.1M.1.counts, add.cell.ids=c("rna.1M.1"),project = "atac.rna.1M.1")  
####1M.2 renames
integarted.atac.1M.2.renames<-merge(integarted.atac.1M.2.counts, add.cell.ids=c("rna.1M.2"),project = "atac.rna.1M.2")  
####2.5M.1 renames
integarted.atac.2.5M.1.renames<-merge(integarted.atac.2.5M.1.counts, add.cell.ids=c("rna.2.5M.1"),project = "atac.rna.2.5M.1")
####2.5M.2 renames
integarted.atac.2.5M.2.renames<-merge(integarted.atac.2.5M.2.counts, add.cell.ids=c("rna.2.5M.2"),project = "atac.rna.2.5M.2")
####2.5M.3 renames
integarted.atac.2.5M.3.renames<-merge(integarted.atac.2.5M.3.counts, add.cell.ids=c("rna.2.5M.3"),project = "atac.rna.2.5M.3")
####2.5M.4 renames
integarted.atac.2.5M.4.renames<-merge(integarted.atac.2.5M.4.counts, add.cell.ids=c("rna.2.5M.4"),project = "atac.rna.2.5M.4")
####3.5M.1 renames
integarted.atac.3.5M.1.renames<-merge(integarted.atac.3.5M.1.counts, add.cell.ids=c("rna.3.5M.1"),project = "atac.rna.3.5M.1")
####3.5M.2 renames
integarted.atac.3.5M.2.renames<-merge(integarted.atac.3.5M.2.counts, add.cell.ids=c("rna.3.5M.2"),project = "atac.rna.3.5M.2")
####4.5M.1 renames
integarted.atac.4.5M.1.renames<-merge(integarted.atac.4.5M.1.counts, add.cell.ids=c("rna.4.5M.1"),project = "atac.rna.4.5M.1")
####4.5M.2 renames
integarted.atac.4.5M.2.renames<-merge(integarted.atac.4.5M.2.counts, add.cell.ids=c("rna.4.5M.2"),project = "atac.rna.4.5M.2")
####RC renames
integarted.atac.RC.renames<-merge(integarted.atac.RC.counts, add.cell.ids=c("rna.RC"),project = "atac.rna.RC")
####WT1 renames
integarted.atac.WT.1.renames<-merge(integarted.atac.WT.1.counts, add.cell.ids=c("rna.WT.1"),project = "atac.rna.WT.1")
####WT2 renames
integarted.atac.WT.2.renames<-merge(integarted.atac.WT.2.counts, add.cell.ids=c("rna.WT.2"),project = "atac.rna.WT.2")
####2W1 renames
integarted.atac.2W.1.renames<-merge(integarted.atac.2W.1.counts, add.cell.ids=c("rna.2W.1"),project = "atac.rna.2W.1")
####2W2 renames
integarted.atac.2W.2.renames<-merge(integarted.atac.2W.2.counts, add.cell.ids=c("rna.2W.2"),project = "atac.rna.2W.2")
}

rm(integarted.atac.1M.1.counts,integarted.atac.1M.2.counts,integarted.atac.2.5M.1.counts,integarted.atac.2.5M.2.counts,integarted.atac.2.5M.3.counts,integarted.atac.2.5M.4.counts,integarted.atac.2W.1.counts,integarted.atac.2W.2.counts,integarted.atac.3.5M.1.counts,integarted.atac.3.5M.2.counts,integarted.atac.4.5M.1.counts,integarted.atac.4.5M.2.counts,integarted.atac.RC.counts,integarted.atac.WT.1.counts,integarted.atac.WT.2.counts)
gc()
rm(integarted.atac.1M.1.renames,integarted.atac.1M.2.renames,integarted.atac.2.5M.1.renames,integarted.atac.2.5M.2.renames,integarted.atac.2.5M.3.renames,integarted.atac.2.5M.4.renames,integarted.atac.2W.1.renames,integarted.atac.2W.2.renames,integarted.atac.3.5M.1.renames,integarted.atac.3.5M.2.renames,integarted.atac.4.5M.1.renames,integarted.atac.4.5M.2.renames,integarted.atac.RC.renames,integarted.atac.WT.1.renames,integarted.atac.WT.2.renames)
gc()

atac.prostate.merge<-merge(integarted.atac.WT.1.renames.merge,y=c(integarted.atac.WT.2.renames.merge,integarted.atac.2W.1.renames.merge,integarted.atac.2W.2.renames.merge,
                                                                  integarted.atac.1M.1.renames.merge,integarted.atac.1M.2.renames.merge,integarted.atac.2.5M.1.renames.merge,
                                                                  integarted.atac.2.5M.2.renames.merge,integarted.atac.2.5M.3.renames.merge,integarted.atac.2.5M.4.renames.merge,
                                                                  integarted.atac.3.5M.1.renames.merge,integarted.atac.3.5M.2.renames.merge,integarted.atac.4.5M.1.renames.merge,
                                                                  integarted.atac.4.5M.2.renames.merge,integarted.atac.RC.renames.merge),project = "ATAC.prostate.combined") 
rm(integarted.atac.1M.1.renames.merge,integarted.atac.1M.2.renames.merge,integarted.atac.2.5M.1.renames.merge,integarted.atac.2.5M.2.renames.merge,integarted.atac.2.5M.3.renames.merge,integarted.atac.2.5M.4.renames.merge,integarted.atac.2W.1.renames.merge,integarted.atac.2W.2.renames.merge,integarted.atac.3.5M.1.renames.merge,integarted.atac.3.5M.2.renames.merge,integarted.atac.4.5M.1.renames.merge,integarted.atac.4.5M.2.renames.merge,integarted.atac.RC.renames.merge,integarted.atac.WT.1.renames.merge,integarted.atac.WT.2.renames.merge)
gc()


















