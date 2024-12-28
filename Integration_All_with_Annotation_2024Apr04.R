# Set WD
setwd("/Volumes/T7/Yaqing")

# Set seed
set.seed(2)

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)
library(SeuratWrappers)
library(reticulate)
library(stringr)

# load data

load("/Volumes/T7/Yaqing/TK_Annotated_V1.Robj")
load("/Volumes/T7/Yaqing/WT_Annotated_V1.Robj")

# merge, JoinLayers, then integrate by Sample (Condition in my case)

All_Annotated_V1 <- merge(TK_Annotated_V1,WT_Annotated_V1)
All_Annotated_V1 <- JoinLayers(All_Annotated_V1)

All_Annotated_V1[["RNA"]] <- split(All_Annotated_V1[["RNA"]], f = All_Annotated_V1$Condition)

All_Annotated_V1 <- NormalizeData(All_Annotated_V1)
All_Annotated_V1 <- FindVariableFeatures(All_Annotated_V1)
All_Annotated_V1 <- ScaleData(All_Annotated_V1)
All_Annotated_V1 <- RunPCA(All_Annotated_V1,npcs=100)

All_Annotated_V1 <- IntegrateLayers(object = All_Annotated_V1, method=CCAIntegration,
                                    orig.reduction="pca",new.reduction="integrated.cca",
                                    verbose = FALSE)

save(All_Annotated_V1,file='All_Annotated_V1.Robj')

# Analysis starts from here......

load('/Volumes/T7/Yaqing/All_Annotated_V1.Robj')

All_Annotated_V1 <- FindNeighbors(All_Annotated_V1,reduction='integrated.cca',dims=1:40)
All_Annotated_V1 <- FindClusters(All_Annotated_V1,resolution=0.4,cluster.name='cca_clusters')
All_Annotated_V1 <- RunUMAP(All_Annotated_V1,reduction='integrated.cca',dims=1:40,reduction.name='umap.cca')


# Generate some plots
png("CellType_DimPlot.png",width = 15, height = 20, units= 'in',res=300)
p1 <- DimPlot(All_Annotated_V1,label=T)
p2 <- DimPlot(All_Annotated_V1,group.by='CellType',label=T)
p3 <- DimPlot(All_Annotated_V1,group.by='CellType',label=T,split.by='SampleType')
p <- plot_grid(p1,p2,p3,nrow=3)
print(p)
dev.off()
# DimPlot(All_Annotated_V1,label=T,split.by='Condition')
# DimPlot(All_Annotated_V1,label=T,group.by='CellType',split.by='SampleType') 

All_Annotated_V1 <- JoinLayers(All_Annotated_V1)
mark_All <- FindAllMarkers(All_Annotated_V1,only.pos=T,min.pct = 0.5,logfc.threshold = 0.25)
mark_All$ratio <- mark_All$pct.1/mark_All$pct.2
mark_All$power <- mark_All$ratio*mark_All$avg_log2FC
View(mark_All)

save(mark_All,file='Yaqing_All_Annotated_V1_markers_2024-04-04.Robj')
write.table(mark_All,file = 'Yaqing_All_Annotated_V1_markers_2024-04-04.csv',sep = ',',row.names = T,col.names = NA)

FeaturePlot(All_Annotated_V1,'Thbs2',label=T) # Fishing in Fibroblasts
FeaturePlot(All_Annotated_V1,'Dcn',label=T) # Cluster 0, 3, Dcn+
FeaturePlot(All_Annotated_V1,'Ccl8',label=T) # Ccl8 Very Specific to 0
FeaturePlot(All_Annotated_V1,'Dkk2',label=T) # Cluster 6, Dkk2+
FeaturePlot(All_Annotated_V1,'Rspo3',label=T)
FeaturePlot(All_Annotated_V1,'Frzb',label=T) # Cluster 2,8 Rspo3+_Frzb+_Fibroblasts
FeaturePlot(All_Annotated_V1,'Top2a',label=T) # 5,7 cycling, 7 Schwann-Cycling, 5 Fibroblasts-Cycling

png("Schwann_Cell_FeaturePlot.png",width = 24, height = 12, units= 'in',res=300) # Cluster 1,4 Schwann-Cells
p1 <- FeaturePlot(All_Annotated_V1,'Sox10',label=T)
p2 <- FeaturePlot(All_Annotated_V1,'Gap43',label=T) # This is cluster 1 specific
p3 <- FeaturePlot(All_Annotated_V1,'Sox2',label=T) 
p4 <- FeaturePlot(All_Annotated_V1,'Ngfr',label=T) 
p5 <- FeaturePlot(All_Annotated_V1,'S100b',label=T) # This is cluster 1 specific
p6 <- FeaturePlot(All_Annotated_V1,'Egr2',label=T) 
p7 <- FeaturePlot(All_Annotated_V1,'Pou3f1',label=T) 
p8 <- FeaturePlot(All_Annotated_V1,'Mbp',label=T) 
p <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2)
print(p)
dev.off()


FeaturePlot(All_Annotated_V1,'Epcam',label=T) # Cluster 11 contains Epcam
FeaturePlot(All_Annotated_V1,'Ptprc',label=T) # Cluster 9 Immune
FeaturePlot(All_Annotated_V1,'Cdh5',label=T) # Cluster 12 Endothelial


FeaturePlot(All_Annotated_V1,'Taco1',label=T) # Cluster 10 Taco1+
FeaturePlot(All_Annotated_V1,'Myh11',label=T) # Cluster 14
FeaturePlot(All_Annotated_V1,'Ebf2',label=T) # part of Cluster 11 
FeaturePlot(All_Annotated_V1,'Ifit3',label=T) # Cluster 13 IFN-Stimulated

DimPlot(All_Annotated_V1)

All_Annotated_V1 <- RenameIdents(All_Annotated_V1,
                                 '0'='Dcn+_Fibroblasts',
                                 '1'='Schwann_Cells',
                                 '2'='Rspo3+_Frzb+_Fibroblasts',
                                 '3'='Dcn+_Fibroblasts',
                                 '4'='Schwann_Cells',
                                 '5'='Fibroblasts_Cycling',
                                 '6'='Dkk2+_Fibroblasts',
                                 '7'='Schwann_Cycling',
                                 '8'='Rspo3+_Frzb+_Fibroblasts',
                                 '9'='Immune',
                                 '10'='Taco1+_Fibroblasts',
                                 '11'='TBD',
                                 '12'='Endothelial',
                                 '13'='IFN_Stimulated_Fibroblasts',
                                 '14'='Myh11+_Fibroblasts')
All_Annotated_V1$Global_CellType <- Idents(All_Annotated_V1)
DimPlot(All_Annotated_V1)

# Subset cluster 11
Cluster_11 <- subset(All_Annotated_V1,idents='TBD')

# Integrate cluster 11
Cluster_11 <- JoinLayers(Cluster_11)

Cluster_11[["RNA"]] <- split(Cluster_11[["RNA"]], f = Cluster_11$Condition)

Cluster_11 <- NormalizeData(Cluster_11)
Cluster_11 <- FindVariableFeatures(Cluster_11)
Cluster_11 <- ScaleData(Cluster_11)
Cluster_11 <- RunPCA(Cluster_11,npcs=100)

# Integrate the 'layers' (condition groups) (Method CCA)
integrated_Cluster_11.cca <- IntegrateLayers(object = Cluster_11, method=CCAIntegration,
                                             orig.reduction="pca",new.reduction="integrated.cca",
                                             verbose = FALSE,k.weight=40) # here add a k.weight=40, because default 100 gives error


integrated_Cluster_11.cca <- FindNeighbors(integrated_Cluster_11.cca,reduction='integrated.cca',dims=1:20)
integrated_Cluster_11.cca <- FindClusters(integrated_Cluster_11.cca,resolution=0.2,cluster.name='cca_clusters')
integrated_Cluster_11.cca <- RunUMAP(integrated_Cluster_11.cca,reduction='integrated.cca',dims=1:20,reduction.name='umap.cca')

# Rename Idents and get some plots for sub cluster 11

Cluster_11_integrated <- JoinLayers(integrated_Cluster_11.cca)
DimPlot(Cluster_11_integrated,label=T)

Cluster_11_integrated <- RenameIdents(Cluster_11_integrated,
                                      '0'='Dcn+_Fibroblasts',
                                      '1'='Adipogenic_Precursor_Cells',
                                      '2'='Epithelial')

Cluster_11_integrated$Global_CellType <- Idents(Cluster_11_integrated)

# sub object
DimPlot(Cluster_11_integrated)
# Cluster_11_integrated$Global_CellType

# barcodes in the sub object Cluster 11
barcodes.to.edit <- colnames(Cluster_11_integrated)

# check idents slot
table(Idents(All_Annotated_V1))


# Convert Global_CellType in whole object to character
All_Global_CellType <- subset(All_Annotated_V1@meta.data, select = 'Global_CellType')
sapply(All_Global_CellType, class)
tmp <- sapply(All_Global_CellType, is.factor)
All_Global_CellType[tmp] <- lapply(All_Global_CellType[tmp], as.character)


# Convert Global_CellType in sub object Cluster 11 to character
Cluster_11_CellType <- subset(Cluster_11_integrated@meta.data, select = 'Global_CellType')
sapply(Cluster_11_CellType, class)
tmp <- sapply(Cluster_11_CellType, is.factor)
Cluster_11_CellType[tmp] <- lapply(Cluster_11_CellType[tmp], as.character)

# Map over
All_Global_CellType[barcodes.to.edit,] <- Cluster_11_CellType

# Check updated Global_CellType
table(All_Global_CellType$Global_CellType)

# Convert character to factor
sapply(All_Global_CellType, class)
tmp <- sapply(All_Global_CellType, is.character)
All_Global_CellType[tmp] <- lapply(All_Global_CellType[tmp], as.factor)

# Rename Global_CellType to Global_CellType_V2
names(All_Global_CellType)[names(All_Global_CellType) == "Global_CellType"] <- "Global_CellType_V2"

# Add Global_CellType_V2 to whole object's meta data
All_Annotated_V1 <- AddMetaData(All_Annotated_V1, All_Global_CellType)

# Get some plots!

DimPlot(All_Annotated_V1,label=T,group.by='Global_CellType')
DimPlot(All_Annotated_V1,label=T,group.by='Global_CellType_V2')

All_Annotated_V2 <- All_Annotated_V1
save(All_Annotated_V2,file='All_Annotated_V2.Robj')


# New Analysis starts from here......

load("/Volumes/T7/Yaqing/All_Annotated_V2.Robj")

All_Annotated_V2 <- FindNeighbors(All_Annotated_V2,reduction='integrated.cca',dims=1:40)
All_Annotated_V2 <- FindClusters(All_Annotated_V2,resolution=0.4,cluster.name='cca_clusters')
All_Annotated_V2 <- RunUMAP(All_Annotated_V2,reduction='integrated.cca',dims=1:40,reduction.name='umap.cca')

png("Global_CellType_DimPlot.png",width = 15, height = 20, units= 'in',res=300)
p1 <- DimPlot(All_Annotated_V2,label=T)
p2 <- DimPlot(All_Annotated_V2,group.by='Global_CellType_V2',label=T)
p3 <- DimPlot(All_Annotated_V2,group.by='Global_CellType_V2',label=T,split.by='SampleType')
p <- plot_grid(p1,p2,p3,nrow=3)
print(p)
dev.off()

All_Annotated_V2 <- SetIdent(All_Annotated_V2,value='Global_CellType_V2')
mark_All <- FindAllMarkers(All_Annotated_V2,only.pos=T,min.pct = 0.5,logfc.threshold = 0.25)
mark_All$ratio <- mark_All$pct.1/mark_All$pct.2
mark_All$power <- mark_All$ratio*mark_All$avg_log2FC
View(mark_All)

save(mark_All,file='Yaqing_Global_Markers_2024-04-04.Robj')



