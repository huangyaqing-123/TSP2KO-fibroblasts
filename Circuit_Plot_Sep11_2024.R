# Use iGraph to get circuit plot:

# First set wd and load our data:

setwd("/Volumes/T7/Yaqing")

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
library(ggpol)
library(poolr)

# Load data
load('/Volumes/T7/Yaqing/fibro.integrated.Robj')
load('/Volumes/T7/Yaqing/All_Annotated_V2.Robj')

# Try "Pdgfb—Pdgfrb" as test Ligand-Receptor mechanism, Schwann-Dcn as edge in the circuit plot
CellType_List <- rownames(table(All_Annotated_V2$Global_CellType_V2))
non_Fibro_CellType_List <- c('Epithelial','Endothelial','Immune')

# This returns a table with node weight (CellType proportion as weight), you can use this to define vertex size
table(All_Annotated_V2$Global_CellType_V2)
sum(table(All_Annotated_V2$Global_CellType_V2))
table(All_Annotated_V2$Global_CellType_V2)/sum(table(All_Annotated_V2$Global_CellType_V2))

feature <- c('Bmp4—Bmpr2')


Data_for_CircuitPlot_1 <- function(feature){
  
  graph_df <- data.frame()
  for(sendingType in CellType_List){
    if(sendingType %in% non_Fibro_CellType_List) next
    sub1 <- subset(fibro.integrated,subset = SendingType == sendingType)
  for(receivingType in CellType_List){
      
      sub2 <- subset(sub1,subset = ReceivingType == receivingType)
      
      if(table(sub2$SampleType)[1] == 0 || table(sub2$SampleType)[2] == 0) next
      
      sub2_TK <- subset(sub2,subset = SampleType == 'Tsp2.Knockout')
      sub2_WT <- subset(sub2,subset = SampleType == 'WildType')
      
      assay <- 'CellToCell'
      
      selected_TK <- rownames(sub2_TK[[assay]]) %in% feature
      if(sum(selected_TK)>0){
        edge.object <- data.frame(feature.value = sub2_TK[[assay]]@layers$data[selected_TK,])
        edge.object$barcode.pair <- rownames(edge.object)
        weight.1 <- mean(edge.object$feature.value)
        weight.2 <- length(edge.object$feature.value)
        TK_edge <- c(sendingType,receivingType,weight.1,weight.2,'TK')
        graph_df <- rbind(graph_df,TK_edge)
      }
      
      selected_WT <- rownames(sub2_WT[[assay]]) %in% feature
      if(sum(selected_WT)>0){
        selected_WT <- rownames(sub2_WT[[assay]]) %in% feature
        edge.object <- data.frame(feature.value = sub2_WT[[assay]]@layers$data[selected_WT,])
        edge.object$barcode.pair <- rownames(edge.object)
        weight.1 <- mean(edge.object$feature.value)
        weight.2 <- length(edge.object$feature.value)
        WT_edge <- c(sendingType,receivingType,weight.1,weight.2,'WT')
        graph_df <- rbind(graph_df,WT_edge)
      }
    }  
    
  }
  
  colnames(graph_df) <- c('SendingType','ReceivingType','mean','count','SampleType')
  
  return(graph_df)
}

# Sample use for the function
# Notice that FOI must be a vector! Even for one mechanism, use c('Pdgfb—Pdgfrb') instead of 'Pdgfb—Pdgfrb'

FOI <- c('Bmp4—Bmpr2')
to_plot_1 <- Data_for_CircuitPlot(c('Pdgfb—Pdgfrb'))

