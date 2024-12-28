setwd("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing")

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


tk1 <- Read10X("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/TK1/raw_feature_bc_matrix")
wt1 <- Read10X("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/WT1/raw_feature_bc_matrix")

tk1.ge <- tk1$`Gene Expression`
wt1.ge <- wt1$`Gene Expression`

tk1.seurat <- CreateSeuratObject(tk1.ge,min.cells = 3,min.features = 200) # Be careful with these thresholds. I know that min.features = 200 is OK here, because this 10X v3 data which has very high UMI
wt1.seurat <- CreateSeuratObject(wt1.ge,min.cells = 3,min.features = 200)


tk1.seurat$Condition <- 'TK1'
wt1.seurat$Condition <- 'WT1'

tk1.seurat$SampleType <- 'Tsp2.Knockout'
wt1.seurat$SampleType <- 'WildType'

tk1.seurat[["percent.mt"]] <- PercentageFeatureSet(tk1.seurat, pattern = "^mt-") # Mouse: "^mt-" | Rat: "^Mt-" | Human: "^MT-"
wt1.seurat[["percent.mt"]] <- PercentageFeatureSet(wt1.seurat, pattern = "^mt-")

tk2 <- Read10X("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/TK2/raw_feature_bc_matrix")
wt2 <- Read10X("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/WT2/raw_feature_bc_matrix")

tk2.seurat <- CreateSeuratObject(tk2,min.cells = 3,min.features = 100) # Be careful with these thresholds. I know that min.features = 200 is OK here, because this 10X v3 data which has very high UMI
wt2.seurat <- CreateSeuratObject(wt2,min.cells = 3,min.features = 100)

tk2.seurat$Condition <- 'TK2'
wt2.seurat$Condition <- 'WT2'
tk2.seurat$SampleType <- 'Tsp2.Knockout'
wt2.seurat$SampleType <- 'WildType'


tk2.seurat[["percent.mt"]] <- PercentageFeatureSet(tk2.seurat, pattern = "^mt-") # Mouse: "^mt-" | Rat: "^Mt-" | Human: "^MT-"
wt2.seurat[["percent.mt"]] <- PercentageFeatureSet(wt2.seurat, pattern = "^mt-")

png('First_Look_QC_TK2.png',width = 10,height = 8,units = 'in',res=300)
p1 <- VlnPlot(tk2.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
p2 <- VlnPlot(tk2.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
print(plot_grid(p1,p2))
dev.off()

png('First_Look_QC_WT2.png',width = 10,height = 8,units = 'in',res=300)
p1 <- VlnPlot(wt2.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
p2 <- VlnPlot(wt2.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
print(plot_grid(p1,p2))
dev.off()


tk3 <- Read10X("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/TK3/raw_feature_bc_matrix")
wt3 <- Read10X("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/WT3/raw_feature_bc_matrix")

tk3.seurat <- CreateSeuratObject(tk3,min.cells = 3,min.features = 100) # Be careful with these thresholds. I know that min.features = 200 is OK here, because this 10X v3 data which has very high UMI
wt3.seurat <- CreateSeuratObject(wt3,min.cells = 3,min.features = 100)

tk3.seurat$Condition <- 'TK3'
wt3.seurat$Condition <- 'WT3'
tk3.seurat$SampleType <- 'Tsp2.Knockout'
wt3.seurat$SampleType <- 'WildType'

tk3.seurat[["percent.mt"]] <- PercentageFeatureSet(tk3.seurat, pattern = "^mt-") # Mouse: "^mt-" | Rat: "^Mt-" | Human: "^MT-"
wt3.seurat[["percent.mt"]] <- PercentageFeatureSet(wt3.seurat, pattern = "^mt-")

png('First_Look_QC_TK3.png',width = 10,height = 8,units = 'in',res=300)
p1 <- VlnPlot(tk3.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
p2 <- VlnPlot(tk3.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
print(plot_grid(p1,p2))
dev.off()

png('First_Look_QC_WT3.png',width = 10,height = 8,units = 'in',res=300)
p1 <- VlnPlot(wt3.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
p2 <- VlnPlot(wt3.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
print(plot_grid(p1,p2))
dev.off()

merge2 <- merge(tk2.seurat,wt2.seurat)
merge3 <- merge(tk3.seurat,wt3.seurat)
merge2_sub <- subset(merge2,subset =
                    nFeature_RNA > 300 &
                    percent.mt < 10 &
                    percent.mt > 0.1)
merge3_sub <- subset(merge3,subset =
                    nFeature_RNA > 300 &
                    percent.mt < 10 &
                    percent.mt > 0.1)

png("First_Look_QC_filtered_TK2_WT2.png",width = 10, height = 8, units= 'in',res=300)
p1 <- VlnPlot(merge2_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
p2 <- VlnPlot(merge2_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T,group.by = 'Condition')
print(plot_grid(p1,p2))
dev.off()

png("First_Look_QC_filtered_TK3_WT3.png",width = 10, height = 8, units= 'in',res=300)
p1 <- VlnPlot(merge3_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T)
p2 <- VlnPlot(merge3_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,log = T,group.by = 'Condition')
print(plot_grid(p1,p2))
dev.off()

TK1_sub_1 <- subset(tk1.seurat,subset = 
                  nFeature_RNA > 1000 &
                  percent.mt < 10 &
                  percent.mt > 0.1)
WT1_sub_1 <- subset(wt1.seurat,subset = 
                      nFeature_RNA > 1000 &
                      percent.mt < 10 &
                      percent.mt > 0.1)


# Filter the data first pass, be careful
TK2_sub_1 <- subset(tk2.seurat,subset = 
                  nFeature_RNA > 300 &
                  percent.mt < 10 &
                  percent.mt > 0.1)
WT2_sub_1 <- subset(wt2.seurat,subset = 
                  nFeature_RNA > 300 &
                  percent.mt < 10 &
                  percent.mt > 0.1)

TK3_sub_1 <- subset(tk3.seurat,subset = 
                      nFeature_RNA > 300 &
                      percent.mt < 10 &
                      percent.mt > 0.1)
WT3_sub_1 <- subset(wt3.seurat,subset = 
                      nFeature_RNA > 300 &
                      percent.mt < 10 &
                      percent.mt > 0.1)

# Run PCA analysis, check the image output to determine dim_num for clustering
PCA_analysis <- function(data,sample_name){
  data <- NormalizeData(data)
  data <- ScaleData(data)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, npcs = 100)
  pdf(file = paste(sample_name,'.PCs.pdf',sep=''), width=10,height=8)
  print(ElbowPlot(data,ndims = 100))
  print(PCHeatmap(data,cells=200,balanced=T,dims=1:9))
  print(PCHeatmap(data,cells=200,balanced=T,dims=10:18))
  print(PCHeatmap(data,cells=200,balanced=T,dims=19:27))
  print(PCHeatmap(data,cells=200,balanced=T,dims=28:36))
  print(PCHeatmap(data,cells=200,balanced=T,dims=37:45))
  print(PCHeatmap(data,cells=200,balanced=T,dims=46:54))
  print(PCHeatmap(data,cells=200,balanced=T,dims=55:63))
  print(PCHeatmap(data,cells=200,balanced=T,dims=64:72))
  print(PCHeatmap(data,cells=200,balanced=T,dims=73:81))
  print(PCHeatmap(data,cells=200,balanced=T,dims=82:90))
  print(PCHeatmap(data,cells=200,balanced=T,dims=91:99))
  dev.off()
  return(data)
}


TK1_sub_1 <- PCA_analysis(TK1_sub_1,'TK1')
WT1_sub_1 <- PCA_analysis(WT1_sub_1,'WT1')
TK2_sub_1 <- PCA_analysis(TK2_sub_1,'TK2')
WT2_sub_1 <- PCA_analysis(WT2_sub_1,'WT2')
TK3_sub_1 <- PCA_analysis(TK3_sub_1,'TK3')
WT3_sub_1 <- PCA_analysis(WT3_sub_1,'WT3')

# Embed and cluster
# dim_num is the number of principal components we want to keep in our further analysis.
# A reasonable dim_num can be estimated by checking the PCs.pdf file, especially the Elbow plot part
# res is the resolution parameter in FindClusters
Embed_and_Cluster <- function(data,data_name,dim_num,res){
  data <- RunUMAP(data, reduction = "pca", dims = 1:dim_num)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:dim_num)
  data <- FindClusters(data, resolution = res)
  pdf(paste(data_name,"_QC_check_after_PCA.pdf",sep=''),width = 10,height = 8)
  print(DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE))
  print(FeaturePlot(data,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)) # Check if there are obvious multi-plet clusters
  print(FeaturePlot(data,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)) # Check which clusters look low-quality/low-information.
  print(VlnPlot(data,c('nFeature_RNA','nCount_RNA','percent.mt')))
  dev.off()
  return(data)
}

# Run Embedding and Clustering for each sample, check QC and clusters

# TK1 (Redo actually, since last time we combined TK1, WT1 together)
TK1_sub_1 <- Embed_and_Cluster(TK1_sub_1,'TK1',20,0.2)
png("QC_check_after_PCA_TK1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK1_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK1.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK1_sub_1, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK1_sub_1,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK1_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()


# WT1 (Redo actually, since last time we combined TK1, WT1 together)
WT1_sub_1 <- Embed_and_Cluster(WT1_sub_1,'WT1',20,0.2)
png("QC_check_after_PCA_WT1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT1_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT1.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT1_sub_1, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_1,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

# Further cleaning for TK1

TK1_sub_2 <- subset(TK1_sub_1,idents = c('1','5'),invert=T)

TK1_sub_2 <- PCA_analysis(TK1_sub_2,'TK1_round1')

TK1_sub_2 <- Embed_and_Cluster(TK1_sub_2,'TK1_round1',20,0.2)

png("QC_check_after_PCA_TK1_round1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK1_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK1_round1.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK1_sub_2, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK1_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK1_sub_2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK1_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK1_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()


# Further cleaning for WT1

WT1_sub_2 <- subset(WT1_sub_1,idents = c('0','1','5'),invert=T)

WT1_sub_2 <- PCA_analysis(WT1_sub_2,'WT1_round1')

WT1_sub_2 <- Embed_and_Cluster(WT1_sub_2,'WT1_round1',20,0.2)

png("QC_check_after_PCA_WT1_round1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT1_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT1_round1.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT1_sub_2, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT1_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT1_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

# Further cleaning for WT1 second pass

WT1_sub_3 <- subset(WT1_sub_2,idents = c('3'),invert=T)

WT1_sub_3 <- PCA_analysis(WT1_sub_3,'WT1_round2')

WT1_sub_3 <- Embed_and_Cluster(WT1_sub_2,'WT1_round2',20,0.2)

png("QC_check_after_PCA_WT1_round2.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT1_sub_3,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT1_round2.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT1_sub_3, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT1_round2.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_3,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT1_round2.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_3,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()


WT1_sub_4 <- subset(WT1_sub_3,idents = c('3'),invert=T)

WT1_sub_4 <- PCA_analysis(WT1_sub_4,'WT1_round2')

WT1_sub_4 <- Embed_and_Cluster(WT1_sub_4,'WT1_round2',20,0.2)

png("QC_check_after_PCA_WT1_round3.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT1_sub_4,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT1_round3.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT1_sub_4, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT1_round3.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_4,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT1_round3.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT1_sub_4,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()


# TK2
TK2_sub_1 <- Embed_and_Cluster(TK2_sub_1,'TK2',20,0.2)
png("QC_check_after_PCA_TK2.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK2_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK2.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK2_sub_1, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK2.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK2_sub_1,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK2.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK2_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

# WT2 

WT2_sub_1 <- Embed_and_Cluster(WT2_sub_1,'WT2',20,0.2)
png("QC_check_after_PCA_WT2.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT2_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT2.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT2_sub_1, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT2.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT2_sub_1,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT2.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT2_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

### Cluster 3 is less informative, cluster 7 and 8 are multiplex?

# TK3

TK3_sub_1 <- Embed_and_Cluster(TK3_sub_1,'TK3',20,0.2)
png("QC_check_after_PCA_TK3.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK3_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK3.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK3_sub_1, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK3.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK3_sub_1,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK3.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK3_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

### Cluster 2 is less informative, cluster 5 and 7 are multiplex?

# WT3
WT3_sub_1 <- Embed_and_Cluster(WT3_sub_1,'WT3',20,0.2)
png("QC_check_after_PCA_WT3.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT3_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT3.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT3_sub_1, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT3.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT3_sub_1,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT3.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT3_sub_1,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

### Cluster 0 is less informative

# Further cleaning for TK2

TK2_sub_2 <- subset(TK2_sub_1,idents = c('0'),invert=T)

TK2_sub_2 <- PCA_analysis(TK2_sub_2,'TK2_round1')

TK2_sub_2 <- Embed_and_Cluster(TK2_sub_2,'TK2_round1',20,0.2)

png("QC_check_after_PCA_TK2_round1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK2_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK2_round1.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK2_sub_2, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK2_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK2_sub_2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK2_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK2_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()



# Further cleaning for WT2
# Throw 2 and 3
WT2_sub_2 <- subset(WT2_sub_1,idents = c('2','3'),invert=T)

WT2_sub_2 <- PCA_analysis(WT2_sub_2,'WT2_round1')

WT2_sub_2 <- Embed_and_Cluster(WT2_sub_2,'WT2_round1',20,0.2)

png("QC_check_after_PCA_WT2_round1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT2_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT2_round1.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT2_sub_2, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT2_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT2_sub_2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT2_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT2_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

# Looks good after first round!


# Further cleaning for TK3

TK3_sub_2 <- subset(TK3_sub_1,idents = c('2'),invert=T)

TK3_sub_2 <- PCA_analysis(TK3_sub_2,'TK3_round1')

TK3_sub_2 <- Embed_and_Cluster(TK3_sub_2,'TK3_round1',20,0.2)

png("QC_check_after_PCA_TK3_round1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK3_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK3_round1.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK3_sub_2, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK3_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK3_sub_2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK3_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK3_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

# Cluster 6 is less informative and splitted in the UMAP.

TK3_sub_3 <- subset(TK3_sub_2,idents = c('6'),invert=T)

TK3_sub_3 <- PCA_analysis(TK3_sub_3,'TK3_round2')

TK3_sub_3 <- Embed_and_Cluster(TK3_sub_3,'TK3_round2',20,0.2)

png("QC_check_after_PCA_TK3_round2.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(TK3_sub_3,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_TK3_round2.png",width=10,height=8,units='in',res=300)
print(DimPlot(TK3_sub_3, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_TK3_round2.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK3_sub_3,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_TK3_round2.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK3_sub_3,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

# Further cleaning for WT3

WT3_sub_2 <- subset(WT3_sub_1,idents = c('0','3'),invert=T)

WT3_sub_2 <- PCA_analysis(WT3_sub_2,'WT3_round1')

WT3_sub_2 <- Embed_and_Cluster(WT3_sub_2,'WT3_round1',20,0.2)

png("QC_check_after_PCA_WT3_round1.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT3_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT3_round1.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT3_sub_2, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT3_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT3_sub_2,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT3_round1.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT3_sub_2,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

## Looks good, seems cluster 7 is not Fibroblasts. Does this matter?

WT3_sub_3 <- subset(WT3_sub_2,idents = c('7'),invert=T)

WT3_sub_3 <- PCA_analysis(WT3_sub_3,'WT3_round2')

WT3_sub_3 <- Embed_and_Cluster(WT3_sub_3,'WT3_round2',20,0.2)

png("QC_check_after_PCA_WT3_round2.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(WT3_sub_3,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_WT3_round2.png",width=10,height=8,units='in',res=300)
print(DimPlot(WT3_sub_3, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_WT3_round2.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT3_sub_3,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_WT3_round2.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT3_sub_3,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

merged_data_1 <- merge(TK1_sub_2,WT1_sub_4)
merged_data_2 <- merge(TK2_sub_2,WT2_sub_2)
merged_data_3 <- merge(TK3_sub_3,WT3_sub_3)

temp <- merge(merged_data_1,merged_data_2)
merged_data <- merge(temp,merged_data_3)

save(merged_data,file = 'merged_data.Robj')
load("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/merged_data.Robj")

merged_data <- PCA_analysis(merged_data,'All_Samples')

merged_data <- Embed_and_Cluster(merged_data,'All_Samples',20,0.2)
png("QC_check_after_PCA_All_Samples.png",width = 10,height = 8,units='in',res=300)
print(VlnPlot(merged_data,c('nFeature_RNA','nCount_RNA','percent.mt')))
dev.off()
png("Cluster_check_after_PCA_All_Samples.png",width=10,height=8,units='in',res=300)
print(DimPlot(merged_data, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
png("Feature_Plot_CellType_All_Samples.png",width=10,height=8,units='in',res=300)
FeaturePlot(merged_data,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)
dev.off()
png("Feature_Plot_All_Samples.png",width=10,height=8,units='in',res=300)
FeaturePlot(merged_data,c('nFeature_RNA','nCount_RNA','percent.mt'),label=T)
dev.off()

merged_data <- JoinLayers(merged_data)
mark <- FindAllMarkers(merged_data,only.pos=T)
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio*mark$avg_log2FC
View(mark)
write.table(mark,file = 'merged_data_Markers.csv',sep = ',',row.names = T,col.names = NA)

png("QC_plots_nFeature_RNA_all.png",width=30,height=30,units='in',res=300)
p1 <- VlnPlot(merged_data,c('nFeature_RNA'),pt.size = 0.01,split.by = 'Condition',group.by = 'seurat_clusters')
p2 <- VlnPlot(merged_data,c('nCount_RNA'),pt.size = 0.01,split.by = 'Condition',group.by = 'seurat_clusters')
p3 <- VlnPlot(merged_data,c('percent.mt'),pt.size = 0.01,split.by = 'Condition',group.by = 'seurat_clusters')
plot_grid(p1,p2,p3,nrow=3,ncol=1)
dev.off()

temp.obj_0 <- subset(merged_data,idents = c('0'))
Idents(temp.obj_0) <- temp.obj_0$SampleType
temp.obj_0 <- JoinLayers(temp.obj_0)
temp.mark_0 <- FindAllMarkers(temp.obj_0,only.pos = T)
temp.mark_0$ratio <- temp.mark_0$pct.1/temp.mark_0$pct.2
temp.mark_0$power <- temp.mark_0$ratio*temp.mark_0$avg_log2FC
View(temp.mark_0)
write.table(temp.mark_0,file = 'Cluster0_Markers.csv',sep = ',',row.names = T,col.names = NA)

temp.obj_4 <- subset(merged_data,idents = c('4'))
Idents(temp.obj_4) <- temp.obj_4$SampleType
temp.obj_4 <- JoinLayers(temp.obj_4)
temp.mark_4 <- FindAllMarkers(temp.obj_4,only.pos = T)
temp.mark_4$ratio <- temp.mark_4$pct.1/temp.mark_4$pct.2
temp.mark_4$power <- temp.mark_4$ratio*temp.mark_4$avg_log2FC
View(temp.mark_4)
write.table(temp.mark_4,file = 'Cluster4_Markers.csv',sep = ',',row.names = T,col.names = NA)

png('test.png',width=10,height=8,units='in',res=300)
DimPlot(merged_data, reduction = "umap", label = TRUE, repel = TRUE,split.by = 'SampleType')
dev.off()
png('test1.png',width=20,height=8,units='in',res=300)
DimPlot(merged_data, reduction = "umap", label = TRUE, repel = TRUE,split.by = 'Condition')
dev.off()

sub_0 <- subset(merged_data,idents = c('0'))
sub_4 <- subset(merged_data,idents = c('4'))

png('Cluster0_QC.png',width = 10,height =8, units='in',res=300)
VlnPlot(sub_0,c('nFeature_RNA','nCount_RNA','percent.mt'))+ 
  NoLegend() + plot_annotation(title ='Cluster 0 QC')&
  theme(plot.title = element_text(size = 20,face = 'bold',hjust = 0.5))
dev.off()
png('Cluster0_QC_by_Condition.png',width=10,height=8,units='in',res=300)
VlnPlot(sub_0,c('nFeature_RNA','nCount_RNA','percent.mt'),group.by = 'Condition')+ 
  NoLegend() + plot_annotation(title ='Cluster 0 QC by Condition')&
  theme(plot.title = element_text(size = 20,face = 'bold',hjust = 0.5))
dev.off()

png('Cluster4_QC.png',width = 10,height =8, units='in',res=300)
VlnPlot(sub_4,c('nFeature_RNA','nCount_RNA','percent.mt'))+ 
  NoLegend() + plot_annotation(title ='Cluster 4 QC')&
  theme(plot.title = element_text(size = 20,face = 'bold',hjust = 0.5))
dev.off()
png('Cluster4_QC_by_Condition.png',width=10,height=8,units='in',res=300)
VlnPlot(sub_4,c('nFeature_RNA','nCount_RNA','percent.mt'),group.by = 'Condition')+ 
  NoLegend() + plot_annotation(title ='Cluster 4 QC by Condition')&
  theme(plot.title = element_text(size = 20,face = 'bold',hjust = 0.5))
dev.off()


cell.frac_0 <- table(Idents(sub_0),sub_0$Condition)
cell.frac_4 <- table(Idents(sub_4),sub_4$Condition)

dev.off()
FeaturePlot(merged_data,c('Ifitm1'),label=T)

DimPlot(merged_data,group.by = 'Condition',shuffle=T)
DimPlot(merged_data,group.by = 'seurat_clusters',shuffle=T)

merged_data <- JoinLayers(merged_data)
mark <- FindMarkers(merged_data,ident.1 = '0',min.pct = 0.25)
