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

#################add ATAC assay to RNA assay
rna.All.combined.integrated.filter[["ATAC"]]<-atac.merge
DefaultAssay(rna.All.combined.integrated.filter) <- "ATAC"
rna.All.combined.integrated.filter <- RunHarmony(
  object = rna.All.combined.integrated.filter,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

#########re-compute the UMAP using corrected LSI embeddings
rna.All.combined.integrated.filter <- RunUMAP(rna.All.combined.integrated.filter, dims = 2:50, reduction = 'harmony',reduction.name = "umap.harmony")
DefaultAssay(rna.All.combined.integrated.filter) <- "ATAC"
rna.All.combined.integrated.filter <- FindMultiModalNeighbors(rna.All.combined.integrated.filter, reduction.list = list("pca", "harmony"), dims.list = list(1:50, 2:50))
rna.All.combined.integrated.filter <- RunUMAP(rna.All.combined.integrated.filter, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

p1 <- DimPlot(rna.All.combined.integrated.filter, reduction = "umap",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(rna.All.combined.integrated.filter, reduction = "umap.harmony",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(rna.All.combined.integrated.filter,reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

CairoPNG(file="umap.all.harmony.dims.png",width=1200,height=480)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()## 

#################################
#Renames
rna.All.combined.integrated.filter <- RenameIdents(rna.All.combined.integrated.filter, '0' = 'Luminal','1' = 'Neuroendocrine','2' = 'Mesenchymal_1','3' = 'Neutrophil')
rna.All.combined.integrated.filter <- RenameIdents(rna.All.combined.integrated.filter ,'4' = 'Macrophage', '5' = 'Basal','6' = 'Endothelial_1','7' = 'Seminal vesicle','8' = 'T cell','9' = 'Mesenchymal_2','10' = 'B cell','11' = 'Endothelial_2','12' = 'Neuron')
My_levels <- c('Luminal','Neuroendocrine','Basal','Seminal vesicle',"Mesenchymal_1","Mesenchymal_2","Endothelial_1","Endothelial_2","Neutrophil","Macrophage","T cell","B cell","Neuron")
Idents(rna.All.combined.integrated.filter) <- factor(Idents(rna.All.combined.integrated.filter), levels= My_levels)
DimPlot(rna.All.combined.integrated.filter, reduction = "umap",label = TRUE)
table(Idents(rna.All.combined.integrated.filter))
table(rna.All.combined.integrated.filter$orig.ident)
write.csv(table(Idents(rna.All.combined.integrated.filter)),"Clusters.csv")

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
p1 <- DimPlot(rna.All.combined.integrated.filter,cols = my36colors, reduction = "umap",  label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(rna.All.combined.integrated.filter, cols = my36colors,reduction = "umap.harmony",  label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(rna.All.combined.integrated.filter,cols = my36colors,reduction = "wnn.umap", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN(integrated)")
CairoPNG(file="umap.final.png",width=1800,height=600)

#######Geneactivity
DefaultAssay(rna.All.combined.integrated.filter) <- 'Geneactivity'
CairoPNG(file="fp.Gene_activities.png",width=1200,height=600)
FeaturePlot(
  object = rna.All.combined.integrated.filter,reduction = "wnn.umap",
  features = c("Krt8","Chga","Krt5","Svs5","Col1a2","Cdh5","Ptprc","Plp1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)
dev.off()

########################Visualization for chromatin accessibility on the locus of representative markers 
DefaultAssay(rna.All.combined.integrated.filter) <- "ATAC"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[1]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/WT_1/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[2]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/WT_2/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[3]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/2W_1/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[4]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/2W_2/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[5]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/1M_1/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[6]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/1M_2/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[7]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_1/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[8]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_2/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[9]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_3/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[10]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/2_5M_4/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[11]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/3_5M_1/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[12]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/3_5M_2/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[13]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/4_5M_1/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[14]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/4_5M_2/atac_fragments.tsv.gz"
rna.All.combined.integrated.filter@assays$ATAC@fragments[[15]]@path<-"E:/lifei/TPPRC/sc/analysis/Aggregation/data/RC/atac_fragments.tsv.gz"

CairoPNG(file="Corplot.Krt8.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Krt8', cols = my36colors,features = 'Krt8', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Krt5.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Krt5', features = 'Krt5', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Chga.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Chga', features = 'Chga', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Svs5.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Svs5', features = 'Svs5', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Cdh5.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Cdh5', features = 'Cdh5', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.S100a9.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='S100a9', features = 'S100a9', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.C1qa.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='C1qa', features = 'C1qa', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Cd28.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Cd28', features = 'Cd28', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Bank1.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Bank1', features = 'Bank1', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Cdh19.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Plp1', features = 'Plp1', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Col1a2.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Col1a2', features = 'Col1a2', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
CairoPNG(file="Corplot.Ptprc.png",width=600,height=480)
CoveragePlot(rna.All.combined.integrated.filter,region ='Ptprc', features = 'Ptprc', assay = 'ATAC', expression.assay = 'RNA', peaks = FALSE)
dev.off()
#######################RNA Scaling

DefaultAssay(rna.All.combined.integrated.filter) <- "RNA"
rna.All.combined.integrated.filter <- ScaleData(rna.All.combined.integrated.filter, verbose = FALSE)

##############################Cell composition
library(ggplot2)
table(Idents(rna.All.combined.integrated.filter))
table(rna.All.combined.integrated.filter$timepoint)
cell.prop<-as.data.frame(prop.table(table(Idents(rna.All.combined.integrated.filter), rna.All.combined.integrated.filter$timepoint)))
head(cell.prop)
colnames(cell.prop)<-c("cluster","timepoint","proportion")
write.csv(cell.prop,"cell.prop.csv")

CairoPNG(file="Cellprop.png",width=600,height=480)
ggplot(cell.prop,aes(timepoint,proportion,fill=cluster),)+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=my36colors)
dev.off()

CairoPNG(file="Cellprop_no_WT.png",width=600,height=480)
ggplot(subset(cell.prop,cell.prop$timepoint!="1_WT"),aes(timepoint,proportion,fill=cluster),)+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
dev.off()
######################Definition for RNA markers
DefaultAssay(rna.All.combined.integrated.filter) <- "integrated"
rna.All.combined.integrated.filter.markers <- FindAllMarkers(rna.All.combined.integrated.filter,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
DefaultAssay(rna.All.combined.integrated.filter) <- "RNA"
rna.All.combined.integrated.filter.markers.rna <- FindAllMarkers(rna.All.combined.integrated.filter,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

######################Extraction of luminal and neuroendocrine cells, and with WT cells exclusion, resulting in 40223 cells
rna.All.combined.integrated.filter.prostate<-subset(rna.All.combined.integrated.filter,(wsnn_res.0.04=="0"|wsnn_res.0.04=="1")&timepoint!="1_WT")
rna.All.combined.integrated.filter.prostate#40223
#####################Low-quality cell removal, resulting in 28754 cells for further analyses 
rna.All.combined.integrated.filter.prostate.filter <- subset(x = rna.All.combined.integrated.filter.prostate,nFeature_RNA>2000)
rna.All.combined.integrated.filter.prostate.filter <- FindMultiModalNeighbors(rna.All.combined.integrated.filter.prostate.filter, reduction.list = list("pca", "harmony"), dims.list = list(1:50, 2:50))
rna.All.combined.integrated.filter.prostate.filter <- RunUMAP(rna.All.combined.integrated.filter.prostate.filter, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
rna.All.combined.integrated.filter.prostate.filter <- FindClusters(rna.All.combined.integrated.filter.prostate.filter, resolution = 0.08,graph.name = "wsnn", algorithm = 1, verbose = FALSE)
CairoPNG(file="Umap_prostate_filter_0.08_1.png",width=600,height=480)
DimPlot(rna.All.combined.integrated.filter.prostate.filter,reduction = "wnn.umap", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")
dev.off()
DimPlot(rna.All.combined.integrated.filter.prostate.filter,reduction = "umap.harmony", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")

CairoPNG(file="umap.prostate.filter.split.png",width=1200,height=300)
DimPlot(rna.All.combined.integrated.filter.prostate.filter, split.by = 'timepoint', reduction = "wnn.umap",pt.size = 0.1) + ggplot2::ggtitle("wnn.umap")
dev.off()## 

CairoPNG(file="fp.Lineage.filter.png",width=1200,height=1000)
FeaturePlot(rna.All.combined.integrated.filter.prostate.filter, features = c("rna_Chga","rna_Ar","rna_Krt8","rna_Col1a1"),reduction = "wnn.umap",cols= c("gray", "red"))
dev.off()
FeaturePlot(rna.All.combined.integrated.filter.prostate.filter, features = c("rna_Ar","rna_Mki67","rna_Chga","rna_Ascl1"),max.cutoff = 2,reduction = "wnn.umap",cols= c("gray", "red"))


#########################renames
table(rna.All.combined.integrated.filter.prostate.filter$timepoint)
#table(rna.All.combined.integrated.filter.prostate.filter$wsnn_res.0.08)
#Idents(rna.All.combined.integrated.filter.prostate.filter) <- factor(rna.All.combined.integrated.filter.prostate.filter$wsnn_res.0.08)

rna.All.combined.integrated.filter.prostate.filter <- RenameIdents(rna.All.combined.integrated.filter.prostate.filter, '3' = '1','0' = '2','2' = '3','1' = '4','4' = '5')
My_levels <- c('1','2','3','4',"5")
Idents(rna.All.combined.integrated.filter.prostate.filter) <- factor(Idents(rna.All.combined.integrated.filter.prostate.filter), levels= My_levels)
rna.All.combined.integrated.filter.prostate.filter$newnames<-factor(Idents(rna.All.combined.integrated.filter.prostate.filter), levels= My_levels)

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
CairoPNG(file="1.png",width=600,height=500)
DimPlot(rna.All.combined.integrated.filter.prostate.filter, reduction = "wnn.umap",cols=my36colors[c(23,25,32,24,34)],label = TRUE)
dev.off()














