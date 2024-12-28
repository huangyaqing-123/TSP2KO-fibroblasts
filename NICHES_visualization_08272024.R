# Set WD
setwd("~/project/scRNA seq analysis 04012024/NICHE_07242024")
options(future.globals.maxSize = 2000 * 1024^2)

# Set seed
set.seed(2)

# Packages
require(Seurat)
require(ggplot2)
require(viridis)
require(RColorBrewer)
require(scales)
require(dplyr)
require(circlize)
require(ComplexHeatmap)
require(cowplot)
require(patchwork)
require(SeuratWrappers)
require(reticulate)
require(stringr)
require(ggpol)
require(poolr)
require(ggrepel)
library(pheatmap)

# Load data
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/fibro.integrated.Robj')
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/signaling.markers.Robj')
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/Markers_Differentiate_TK.Robj')
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/Markers_Differentiate_WT.Robj')

# Heatmap
# Calculate CelltoCell average expression
Idents(fibro.integrated) = fibro.integrated@meta.data$Condition
meanmat = AverageExpression(object = fibro.integrated, assays = "CellToCell")
meanmat = meanmat$CellToCell

# Identify pairs
Top_WT = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 7,]
Top_10_WT = Top_WT[order(Top_WT$`Avg avg_log2FC`, decreasing = TRUE),]$gene[1:10]
Top_TK = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 7,]
Top_10_TK = Top_TK[order(Top_TK$`Avg avg_log2FC`, decreasing = TRUE),]$gene[1:10]

# Draw heatmap
marker_genes = c(Top_10_TK, Top_10_WT)
drawmat = t(as.matrix(meanmat[match(marker_genes,rownames(meanmat)),]))
drawmat = apply(drawmat,2,scale)
rownames(drawmat) = colnames(meanmat)
p1<- pheatmap(drawmat, cluster_cols=F, 
         cluster_rows=F, 
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100),
         border_color = NA)


RidgePlot(fibro.integrated, features =  "Ank3", group.by = 'SampleType')

# Cluster featurePlot
DimPlot(fibro.integrated, reduction = 'umap.rpca', label = T) + NoAxes() + NoLegend()
FeaturePlot(fibro.integrated,reduction='umap.rpca', "Cdh1—Egfr", order = T) + NoAxes()
p1<-FeaturePlot(fibro.integrated,reduction='umap.rpca', "Fgf5—Fgfr1", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p2<-FeaturePlot(fibro.integrated,reduction='umap.rpca', "Wnt4—Fzd2", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p2.5<-FeaturePlot(fibro.integrated,reduction='umap.rpca', "Fgf10—Fgfr1", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p3<-FeaturePlot(fibro.integrated,reduction='umap.rpca', "Ereg—Erbb3", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p4 <- FeaturePlot(fibro.integrated,reduction='umap.rpca', "Ngf—Ngfr", order = F) + NoAxes()+
theme(text=element_text(family="Arial", face="bold"))

p5 <- FeaturePlot(fibro.integrated,reduction='umap.rpca', "Thbs1—Cd36", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p6 <- FeaturePlot(fibro.integrated,reduction='umap.rpca', "Efnb2—Pecam1", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p7 <- FeaturePlot(fibro.integrated,reduction='umap.rpca', "Pf4—Fgfr2", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p8 <- FeaturePlot(fibro.integrated,reduction='umap.rpca', "Bmp6—Bmpr2", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))
p9 <- FeaturePlot(fibro.integrated,reduction='umap.rpca', "Cdh1—Egfr", order = F) + NoAxes()+
  theme(text=element_text(family="Arial", face="bold"))

print((p1|p2|p2.5|p3|p4)/(p5|p6|p7|p8|p9))

# Top5 genes for each one
marker_1 <- New_marker[New_marker$cluster == "SCP - Fb", ]
marker_1 <- marker_1[order(marker_1$avg_log2FC, decreasing = T),][1:5,]
#

# Heatmap
# source function 
source("~/project/scRNA seq analysis 04012024/NICHE_07242024/CustomHeatmap.R")
DoHeatmap(fibro.integrated, features = New_marker, group.by = "New_Idents")
CustomHeatmap(fibro.integrated,
                          data.type = 'CellToCell',
                          primary = 'New_Ident' ,
                          secondary = 'New_Ident' ,
                          tertiary = 'New_Ident' ,
                          ### #quarternary = 'orig.ident' ,
                          primary.cols = NULL,
                          secondary.cols = NULL, # Need to be a named list of colors
                          tertiary.cols = NULL,
                          ### #quarternary.cols = NULL,
                          features = New_marker,
                          labels = c('Signalign Types','Signalign Types','Signalign Types'),
                          selected.row.anotations=NULL,
                          selected.label.size = 10,
                          use.scale.data = T,
                          range.frac = 1)

##category
source("~/project/scRNA seq analysis 04012024/NICHE_07242024/RuleSetFunction.R")
source("~/project/scRNA seq analysis 04012024/NICHE_07242024/DefineCategories.R")

TK_marker <- Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg avg_log2FC`>=1,]$gene
WT_marker <- Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg avg_log2FC`>=1,]$gene

DefineCategories(TK_marker, # A vector of mechanisms, consisting of a sending part and receiving part
                             rule.set.function = rule.set.function, # A function containing categorization rules
                             check.for.consistency = TRUE # Whether or not to check if rule.set is self-consistent
) #not complete

## Subset signalings
PDGF <- subset(fibro.integrated, features = c("Pdgfb—Pdgfrb","Pdgfb—Lrp1","Pdgfb—Itgav"))
PDGF_lrp <- subset(fibro.integrated, features = "Pdgfb—Lrp1")
BMP4 <- subset(fibro.integrated, features = c("Bmp4—Acvr1", "Bmp4—Bmpr1a", "Bmp4—Bmpr2"))


## Cell proportion
# Compute
data_pdgf_sending <- prop.table(table(PDGF$SampleType,PDGF$SendingType),1)*100
data_pdgf_receiving <- prop.table(table(PDGF$SampleType,PDGF$ReceivingType),1)*100
# Define total number of cells, to label x-axis bins
data0 <- as.data.frame(table(PDGF$SampleType))
# Organize for plotting
to.plot.sending <- as.data.frame(data_pdgf_sending)
# Make plot
p_PDGF_s <- ggplot(to.plot.sending, 
       aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = colors.use, name="") + 
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("") + 
  ylab("") +
  labs(title = 'Sending Cell Type') +
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size= 16,color='black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  geom_text(aes(x = Var1, y = 1, label = Freq, fill = NULL), color='black',
            data=data0, position = position_dodge(width = 0.5),vjust=-0.3,size=3)
p_PDGF_s 

# -----RidgePlot for Signaling-----
# Rename idents
Idents(All_Annotated_V2) <- "Global_CellType_V2"
All_Annotated_V2 <- RenameIdents(All_Annotated_V2,
                                 'Dcn+_Fibroblasts'='Dcn+_Fb',
                                 'Dkk2+_Fibroblasts' = 'Dkk2+_Fb', 
                                 'Rspo3+_Frzb+_Fibroblasts' = 'Rspo3+_Fb',
                                 'Fibroblasts_Cycling' = 'Fb_Cycling',
                                 'Schwann_Cycling' = 'Sc_Cycling',
                                 'Schwann_Cells' = 'Sc_Precursors',
                                 'IFN_Stimulated_Fibroblasts' = 'IFN_Stimulated_Fb',
                                 'Adipogenic_Precursor_Cells' = 'Ebf2+_Fb',
                                 'Taco1+_Fibroblasts' = 'Taco1+_Fb',
                                 'Myh11+_Fibroblasts' = 'VSMC',
                                 'Immune'='Immune',
                                 'Epithelial' = 'Epithelial',
                                 'Endothelial'='Endothelial')
All_Annotated_V2$Global_CellType_V3 <- Idents(All_Annotated_V2)

# Color 
colors.use <- c('#59C9A5','#F18F01',
                '#E072A4','#5FA8D3',
                '#0267C1', '#B3001B','#A034F0','#C490D1', '#7C6A0A','#053B06',
                'green','#53599A','#F0CF65')
names(colors.use) <- unique(All_Annotated_V2$Global_CellType_V3)

# RidgePlot for PDGFb
RidgePlot(All_Annotated_V2, "Pdgfb", cols = colors.use, split.by = "Sample Type")
VlnPlot(All_Annotated_V2, "Pdgfb", cols = colors.use, pt.size = 0, split.by = "SampleType") + 
  NoLegend()
colors.use.2 <- c("#ae282c","#2066a8")
names(colors.use.2) <- unique(All_Annotated_V2$SampleType)