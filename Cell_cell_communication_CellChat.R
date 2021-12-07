setwd("E:/lifei/TPPRC/sc/analysis/Aggregation/outs/Linux/Integrated/Re-analysis-all/")
{library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(EnsDb.Mmusculus.v79)
  library(Signac)
  library(S4Vectors)
  library(patchwork)
  library(harmony)
  library(Cairo)
  library(slingshot)
  library(TSCAN)
  library(scales)
  library(viridisLite)
  library(DDRTree)
  library(tradeSeq)
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(svglite)}
load("E:/lifei/TPPRC/sc/analysis/Aggregation/outs/Linux/Integrated/Re-analysis-all/rna.All.combined.integrated.filter.20210519.Rdata")
table(rna.All.combined.integrated.filter@meta.data$timepoint)
head(Idents(rna.All.combined.integrated.filter))
rna.All.combined.integrated.filter$cellname<-Idents(rna.All.combined.integrated.filter)
data.input.tpprc <- GetAssayData(rna.All.combined.integrated.filter, assay = "RNA", slot = "data") 
dim(data.input.tpprc)
meta.tpprc <-rna.All.combined.integrated.filter@meta.data
cellchat.tpprc <- createCellChat(object = data.input.tpprc,meta=meta.tpprc,group.by ="cellname")
summary(cellchat.tpprc)
levels(cellchat.tpprc@idents)## show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat.tpprc@idents))## number of cells in each cell group
####
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search =c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))# use Secreted Signaling for cell-cell communication analysis
table(CellChatDB.use$interaction$annotation)
showDatabaseCategory(CellChatDB.use)
#CellChatDB.use.Secreted.Signaling <- subsetDB(CellChatDB, search =c("Secreted Signaling"))# use Secreted Signaling for cell-cell communication analysis
#dim(CellChatDB.use.Secreted.Signaling$interaction)
dim(CellChatDB.use$interaction)
cellchat.tpprc@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

cellchat.tpprc <- subsetData(cellchat.tpprc) # subset the expression data of signaling genes for saving computation cost
  #future::plan("multiprocess", workers = 4) # do parallel  
  cellchat.tpprc <- identifyOverExpressedGenes(cellchat.tpprc)
  cellchat.tpprc <- identifyOverExpressedInteractions(cellchat.tpprc)
  cellchat.tpprc <- projectData(cellchat.tpprc, PPI.mouse)  
  
  cellchat.tpprc <- computeCommunProb(cellchat.tpprc)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat.tpprc <- filterCommunication(cellchat.tpprc, min.cells = 10)
  cellchat.tpprc <- computeCommunProbPathway(cellchat.tpprc)
  cellchat.tpprc <- aggregateNet(cellchat.tpprc)
  groupSize <- as.numeric(table(cellchat.tpprc@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat.tpprc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat.tpprc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  mat <- cellchat.tpprc@net$weight
  par(mfrow = c(3,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  
  ####
  # Compute the network centrality scores
  cellchat.tpprc <- netAnalysis_computeCentrality(cellchat.tpprc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat.tpprc)
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest
  gg2 <- netAnalysis_signalingRole_scatter(cellchat.tpprc, signaling = c("CXCL", "CCL"))
  #> Signaling role analysis on the cell-cell communication network from user's input
  gg1 + gg2
  
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat.tpprc, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat.tpprc, pattern = "incoming")
  ht1 + ht2
##########Systems analysis of cell-cell communication network
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
library(NMF)
library(ggalluvial)
selectK(cellchat.tpprc, pattern = "outgoing")
nPatterns = 3
cellchat.tpprc <- identifyCommunicationPatterns(cellchat.tpprc, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat.tpprc, pattern = "outgoing")
netAnalysis_dot(cellchat.tpprc, pattern = "outgoing")

selectK(cellchat.tpprc, pattern = "incoming")
nPatterns = 4
cellchat.tpprc <- identifyCommunicationPatterns(cellchat.tpprc, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat.tpprc, pattern = "incoming")
netAnalysis_dot(cellchat.tpprc, pattern = "incoming")
################Identify signaling groups based on their functional similarity
save(cellchat.tpprc, file = "cellchat.tpprc.RData")  

setwd("E:/lifei/TPPRC/sc/analysis/Aggregation/outs/Linux/Integrated/Re-analysis-all/")
load("cellchat.tpprc.RData")
library(CellChat)

groupSize <- as.numeric(table(cellchat.tpprc@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.tpprc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.tpprc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.tpprc@net$weight
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pheatmap::pheatmap(mat,cluster_rows = FALSE,cluster_cols = FALSE)
mat_dot<-c(as.numeric(mat[,1]),as.numeric(mat[,2]),as.numeric(mat[,3]),as.numeric(mat[,4]),as.numeric(mat[,5]),as.numeric(mat[,6]),as.numeric(mat[,7]),as.numeric(mat[,8]),as.numeric(mat[,9]),as.numeric(mat[,10]),as.numeric(mat[,11]),as.numeric(mat[,12]),as.numeric(mat[,13]))
mat_dot<-data.frame(mat_dot,rep(colnames(mat),13))
mat_dot$incoming<-c(rep("Luminal",13),rep("Neuroendocrine",13),rep("Basal",13),rep("Seminal vesicle",13),rep("Mesenchymal_1",13),rep("Mesenchymal_2",13),rep("Endothelial_1",13),rep("Endothelial_2",13),rep("Neutrophil",13),rep("Macrophage",13),rep("T cell",13),rep("B cell",13),rep("Neuron",13))
colnames(mat_dot)<-c("Communication_potential","Outgoing","Incoming")
colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
mat_dot$Incoming<-factor(mat_dot$Incoming,levels=unique(mat_dot$Incoming))
mat_dot$Outgoing<-factor(mat_dot$Outgoing,levels=(unique(mat_dot$Incoming)))

ggplot(data = mat_dot,mapping = aes(x=Incoming,y = Outgoing))+
  geom_point(aes(size = Communication_potential,color=Incoming))+scale_color_manual(values=(colorSpace[1:13]))+
  theme_bw()+scale_size(range = c(1, 20))+
  labs(title = 'Communication strength across cell types',
       x = 'Incoming (Recepting signals)',
       y = 'Outgoing (Sending signals)')

ggplot(data = mat_dot,mapping = aes(x=Incoming,y = Outgoing))+
  geom_point(aes(color=Communication_potential,size = Communication_potential))+ scale_colour_viridis_c(option = "viridis")+
  theme_bw()+scale_size(range = c(1, 20))+
  labs(title = 'Communication strength across cell types',
       x = 'Incoming (Recepting signals)',
       y = 'Outgoing (Sending signals)')

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.tpprc, pattern = "outgoing",width = 15, height = 18)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.tpprc, pattern = "incoming",width = 15, height = 18)
ht1 + ht2

library(NMF)
library(ggalluvial)
selectK(cellchat.tpprc, pattern = "incoming")
nPatterns = 6
cellchat.tpprc <- identifyCommunicationPatterns(cellchat.tpprc, pattern = "incoming", k = nPatterns,width = 15, height = 18)
# river plot
netAnalysis_river(cellchat.tpprc, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
p1<-netAnalysis_dot(cellchat.tpprc, pattern = "incoming",dot.size = c(1, 6))
p1
p1.data<-p1$data
p1.data.NE<-p1.data[which(p1.data$CellGroup=="Neuroendocrine"),]
write.csv(p1.data.NE,"p1.data.NE.csv")

############################Visualization on neuroendocrine-specific incoming communications
library(ggpubr)
library(RColorBrewer)
p1.data.NE.1<-as.data.frame(p1.data.NE[which(p1.data.NE$Contribution>0),])
row.names(p1.data.NE.1)<-p1.data.NE.1$Signaling
ggdotchart(p1.data.NE.1, x = "Signaling", y = "Contribution",
           color = "Contribution",
           sorting = "descending",                        
           add = "segments",                             
           xlab="Signalling pathway", 
           rotate = TRUE,
           dot.size = 6 #,
           #label=deg_number$deg_number
           )

p1.data.NE.1.specific<-p1.data.NE.1[c("GDF","ACTIVIN","DESMOSOME","NEGR","MPZ","NRXN","KIT"),]
ggdotchart(p1.data.NE.1.specific, x = "Signaling", y = "Contribution",
           color = "Contribution",
           sorting = "descending",                        
           add = "segments",                             
           xlab="Signalling pathway", 
           rotate = TRUE,
           dot.size = 6 #,
           #label=deg_number$deg_number
)

selectK(cellchat.tpprc, pattern = "outgoing")
nPatterns = 6
cellchat.tpprc <- identifyCommunicationPatterns(cellchat.tpprc, pattern = "outgoing", k = nPatterns,width = 15, height = 18)
# river plot
netAnalysis_river(cellchat.tpprc, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat.tpprc, pattern = "outgoing",dot.size = c(1, 6))

selected.pathways.incoming<-c("NCAM","SEMA6","NRXN","KIT","NEGR","MPZ","L1CAM","ACTIVIN","GDF") 
selected.pathways.outgoing<-c("NCAM","SPP1","PTN","NRXN","NEGR","CALCR","L1CAM","NT")
pathways.show <- "SEMA6"
group.cellType <- c(names(table(cellchat.tpprc@idents)))
names(group.cellType)<-levels(cellchat.tpprc@idents)
netVisual_chord_cell(cellchat.tpprc, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

par(mfrow=c(1,1))
netVisual_aggregate(cellchat.tpprc, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat.tpprc, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat.tpprc, signaling = pathways.show)
plotGeneExpression(cellchat.tpprc, signaling = "KIT")

save(cellchat.tpprc,file="cellchat.tpprc.RData")

