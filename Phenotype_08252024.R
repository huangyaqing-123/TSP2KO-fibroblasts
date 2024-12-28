# Set working directory 
setwd("~/project/scRNA seq analysis 04012024")

# Package required 
require(Seurat)
require(ggplot2)

#load data
load("~/project/scRNA seq analysis 04012024/All_Annotated_V2.Robj")
load("~/project/scRNA seq analysis 04012024/Yaqing_Global_Markers_2024-04-04.Robj")

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
DimPlot(All_Annotated_V2, cols = colors.use, split.by = "SampleType", label = T, repel = T)+NoLegend()+
  xlab("UAMP-1")+
  ylab("UAMP-2")
DimPlot(All_Annotated_V2, cols = colors.use)+
  xlab("UAMP-1")+
  ylab("UAMP-2")


#Color codes
#color codes for 13 clusters
colors.use <- c('#59C9A5','#F18F01',
                '#E072A4','#5FA8D3',
                '#0267C1', '#B3001B','#A034F0','#C490D1', '#7C6A0A','#053B06',
                 'green','#53599A','#F0CF65')
names(colors.use) <- unique(All_Annotated_V2$Global_CellType_V3)



science.color.codes <- 
heatmap.scale.color <- c("#2066a8",  # Dark Blue
                         "#8ec1da",  # Medium Blue
                         "#cde1ec",  # Light Blue
                         "#ededed",  # Grey
                         "#f6d6c2",  # Light Red
                         "#d47264",  # Medium Red
                         "#ae282c"  )
genotype.color <- c("#2066a8",  # Dark Blue
                         "#8ec1da",  # Medium Blue
                         "#cde1ec",  # Light Blue
                         "#f6d6c2",  # Light Red
                         "#d47264",  # Medium Red
                         "#ae282c"  )
names(genotype.color) <- unique(All_Annotated_V2$Condition)

#Idents for different color
names(color.use.4) <- unique(All_Annotated_V2$CellClass)

#Set Idents
Idents(All_Annotated_V2) <- "Global_CellType_V2"
DimPlot(All_Annotated_V2, shuffle = T)
DimPlot(All_Annotated_V1, split.by = "SampleType", shuffle = T)


#Class category
p1 <- FeaturePlot(All_Annotated_V2,'Epcam',order=T) # Check Epithelial
p2 <- FeaturePlot(All_Annotated_V2,'Ptprc',order=T) # Check Immune
p3 <- FeaturePlot(All_Annotated_V2,'Col1a1',order=T) # Check Mesenchyme
p4 <- FeaturePlot(All_Annotated_V2,'Cdh5',order=T) # Check Endothelial
FeaturePlot(All_Annotated_V2,'Sox10',label=T,order=T) # Check Glial
class_marker <- ((p1|p2)/(p3|p4))
class_marker

#Cluster Marker list
Schwann_marker <- c("Ngfr", "Sox10", "Sox2", "Foxd3", "Erbb3", "S100b", "Pax3", 
                    "Gap43", "Plp1", "Nrcam", "Cdh19")
Dcn_marker <- c("Dcn", "Thy1", "Ccl7", "Fmod", "Prelp", "Aspn", "Pdpn", "Adamts2", "Loxl1")
Dkk_marker <- c("Dkk2", "Ptgs2",  "Ptgs1", "Ereg", "Plpp2")
AP_marker <- c("Ebf2", "Dlk1", "Vit", "Igfbp4", "Plscr2")
IFN_marker <- c("Ifit3", "Oasl2", "Ifit1", "Ifit3b", "Isg15")
Rspo_marker <- c("Rspo3", "Frzb", "Lef1", "Ptprk")
Cycling_marker <- c("Top2a", "Ccna2", "Ccnb1", "Cdk1", "Ccne1",  "Mcm6", "Kif11", "Kif4")
Myh11_marker <- c("Myh11", "Sox5", "Flt1", "Myom1", "Ankrd1")
Taco_marker <- c("Taco1", "Adamts12")
total_marker <- c(Schwann_marker,Dcn_marker,Dkk_marker, AP_marker, IFN_marker,
                  Rspo_marker, Cycling_marker, Myh11_marker, Taco_marker)


#Cell cluster plot
DimPlot(All_Annotated_V2, group.by = "CellClass")
DimPlot(All_Annotated_V2, group.by = "Global_CellType_V2", cols = colors.use)
DimPlot(All_Annotated_V2, split.by = "SampleType", label = T, shuffle = T)
DimPlot(All_Annotated_V2, split.by = "SampleType", group.by = "CellType")

#Subset fibroblast population
Fibro <- subset(All_Annotated_V2, subset = CellClass == "Fibroblasts")
fibro <- subset(Fibro, idents = c("Immune", "Endothelial", "Epithelial"), invert = T)


#Cluster identification plots
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(All_Annotated_V2, features = "Dcn", ncol = 2)
RidgePlot(fibro, features = "Bmp4")

# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(All_Annotated_V2, "Sox10")
VlnPlot(fibro, "Sox10", split.by = "SampleType")
VlnPlot(fibro, "Sox2", split.by = "SampleType")W
VlnPlot(fibro, "S100b", split.by = "SampleType")
VlnPlot(fibro, "Tgfb3", split.by = "SampleType")
VlnPlot(fibro, "Wnt4", split.by = "SampleType")
VlnPlot(fibro, "Pdgfb", split.by = "SampleType")
VlnPlot(fibro, "Lrp1", split.by = "SampleType")
VlnPlot(fibro, "Bmp4", split.by = "SampleType")
feature <- "Timp1"
p6 <- FeaturePlot(All_Annotated_V2, features = feature, split.by = "SampleType")
p7 <- VlnPlot(All_Annotated_V2, features = feature, split.by = "SampleType") +
  theme(axis.title.x = element_blank())
p8 <- (p6/p7)
print(p8)

# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(All_Annotated_V2, features = "Dcn")
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(All_Annotated_V2, 
        features = c("Col1a1","Dcn", "Thy1",
                     "Dkk2", "Bmpr1b", "Ptgs2",
                     "Rspo3","En1","Wnt5a",
                     "Top2a",
                     "Sox10","Ngfr","Erbb3",
                     "Ifit3", "Ifit1", "Oasl2",
                     "Ebf2", "Igfbp4", "Vit",
                     "Taco1",
                     "Myh11", "Myom1", "Cnn1",
                     "Ptprc", "Epcam", "Cdh5")) + 
  RotatedAxis()
# Single cell heatmap of feature expression
gradient_col <- scale_fill_gradientn(colors = heatmap.scale.color)
 
DoHeatmap(subset(fibro, downsample = 200), 
          features = total_marker, 
          size = 4) + NoLegend()



##Proportion analysis 
# Compute
data2 <- prop.table(table(All_Annotated_V2$Condition,All_Annotated_V2$Global_CellType_V2),1)*100
rowSums(data2) # should all be 100
data2
# Define total number of cells, to label x-axis bins
data0 <- as.data.frame(table(All_Annotated_V2$Condition))
# Organize for plotting
to.plot <- as.data.frame(data2)
ggplot(to.plot, 
       aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = colors.use,name="") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = ' ') + 
  xlab("") + 
  ylab("% of Individual Sample") +
  labs(title = 'Composition of Sample') + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.text.x=element_text(size=16,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='black',
            data=data0,position = position_dodge(width = 0.5),vjust=-0.3,size=4)+
  scale_color_manual(values = "red")

##Proportion analysis(Celltype/genotype)
data0.5 <- prop.table(table(All_Annotated_V2$SampleType,All_Annotated_V2$Global_CellType_V3),1)*100
data1.5 <- as.data.frame(table(All_Annotated_V2$SampleType))
to.plot0.5 <- as.data.frame(data0.5)
ggplot(to.plot0.5, 
       aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = colors.use,name="") + 
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("") + 
  ylab("% of Genotype") +
  labs(title = 'Cell Type per Group') + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size = 12),
        axis.text.x=element_text(size=14,color='black'),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        axis.line = element_line(size = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='black',
            data=data1.5,position = position_dodge(width = 0.5),vjust=-0.3,size=4)+
  scale_color_manual(values = "red")


##Proportion analysis (sample/celltype)
data3 <- prop.table(table(fibro$Global_CellType_V2,fibro$Condition),1)*100
rowSums(data3)
data3 <- na.omit(data3)
data3
# Define total number of cells, to label x-axis bins
data4 <- as.data.frame(table(fibro$Global_CellType_V2))
data4[data4==0] <- NA
data4 <- na.omit(data4)
to.plot.2 <- as.data.frame(data3)
ggplot(to.plot.2, 
       aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = genotype.color, name="") + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = '') + 
  xlab("") + 
  ylab("% of Cell Types") +
  labs(title = 'Composition of Individual Cell Type') + 
  theme(plot.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.text.x=element_text(size=16,color='black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='black',
            data=data4, position = position_dodge(width = 0.5),vjust=-0.3,size=4)

## Sample/fibroblast cell types
# Compute
data3.5 <- prop.table(table(fibro$Global_CellType_V2,fibro$Condition),1)*100
# Clean data3.5
data3.5 <- data3.5[-(4:5),]
data3.5 <- data3.5[-(6:6),]
# Define total number of cells, to label x-axis bins
data4.5 <- as.data.frame(table(fibro$Global_CellType_V2))
data4.5[data4.5==0] <- NA
data4.5<-data4.5[complete.cases(data4.5),]
# Organize for plotting
to.plot.3.5 <- as.data.frame(data3.5)
ggplot(to.plot.3.5, 
       aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = genotype.color,name="") + 
  scale_y_continuous(labels = scales::percent_format()) +
  xlab("") + 
  ylab("% of CellType") +
  labs(title = 'Sample per CellType') + 
  theme(plot.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.text.x=element_text(size=14,color='black', angle = 45, hjust=1),
        axis.text.y = element_text(size = 14,color='black'),
        axis.title.x = element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        axis.line = element_line(size = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(aes(x = Var1,y = 1,label = Freq,fill = NULL),color='black',
            data=data4.5,position = position_dodge(width = 0.5),vjust=-0.2,size=4)+
  scale_color_manual(values = "red")


## for single feature, draw feature plot + voilin plot splited by genotype
FeaturePlot(All_Annotated_V2, "Sox10", split.by = "SampleType")
FeaturePlot(All_Annotated_V2, "Lef1", split.by = "SampleType", order = T)
p1<-FeaturePlot(All_Annotated_V2, "Pdgfra", order = T)
p2<-FeaturePlot(All_Annotated_V2, "S100a4", order = T)
p3<-FeaturePlot(All_Annotated_V2, "Dlk1", order = T)
p4<-FeaturePlot(All_Annotated_V2, "Ly6a", order = T)
((p1|p2)/(p3|p4))
FeaturePlot(All_Annotated_V2, "En1", split.by = "SampleType", order = T)


## For cluster markers
FeaturePlot(All_Annotated_V2, "Lef1", order = T) # Wnt-related 
FeaturePlot(All_Annotated_V2, "Axin2", order = T) # origin check 
FeaturePlot(All_Annotated_V2, "Dkk1", order = T) # almost no - check wnt 
FeaturePlot(All_Annotated_V2, "Fgf10", order = T)
FeaturePlot(All_Annotated_V2, "Nppc", order = T) # almost no - check muscle origin
FeaturePlot(All_Annotated_V2, "Adrp", order = T) # no - lipofibroblasts
FeaturePlot(All_Annotated_V2, "Klf4") #transcription facotr for AP
FeaturePlot(All_Annotated_V2, "Pparg") # TFs for AP not concentrated
FeaturePlot(All_Annotated_V2, "Cd24") # TFs for AP not concentrated

FeaturePlot(All_Annotated_V2, "Igfbp4", order = T) # expressed in adipoctyes
FeaturePlot(All_Annotated_V2, "Myod1", order = T) # muscle cell marker - not found 
FeaturePlot(All_Annotated_V2, "Cspg4") # myh11 high expression 
FeaturePlot(All_Annotated_V2, "Cnn1")
FeaturePlot(All_Annotated_V2, "Nr1c3")

#----Isolated SCP----
Sc_precursors <- subset(All_Annotated_V2, subset = Global_CellType_V3 == "Sc_Precursors")
Idents(Sc_precursors) <- "SampleType"
colors.use.2 <- c("#ae282c","#2066a8")
names(colors.use.2) <- unique(All_Annotated_V2$SampleType)
VlnPlot(Sc_precursors, "Sox10", pt.size = 0, cols = colors.use.2)
RidgePlot(Sc_precursors, "Sox10", ncol = 2)

# -----RidgePlot for Signaling




