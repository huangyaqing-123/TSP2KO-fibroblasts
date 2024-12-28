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
library(patchwork)
library(SeuratWrappers)
library(reticulate)
library(stringr)


##### SAM EDITS###
# 1. load data
load("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/merged_data.Robj")
merged_data
# 2. split into two objects
merged_data.list <- SplitObject(merged_data,split.by = "SampleType")
names(merged_data.list)
TK <- merged_data.list$Tsp2.Knockout
WT <- merged_data.list$WildType
TK <- JoinLayers(TK)
WT <- JoinLayers(WT)

# 3. split each individual object into three samples (by Condition)
TK[["RNA"]] <- split(TK[["RNA"]], f = TK$Condition)
WT[["RNA"]] <- split(WT[["RNA"]], f = WT$Condition)

# 4. embed each, via individual 'layers'
TK <- NormalizeData(TK)
TK <- FindVariableFeatures(TK)
TK <- ScaleData(TK)
TK <- RunPCA(TK)
WT <- NormalizeData(WT)
WT <- FindVariableFeatures(WT)
WT <- ScaleData(WT)
WT <- RunPCA(WT)

# 4. integrate the 'layers' (condition groups) for each SampleType (Method CCA)
integrated_TK.cca <- IntegrateLayers(object = TK, method=CCAIntegration,
                                 orig.reduction="pca",new.reduction="integrated.cca",
                                 verbose = FALSE)

integrated_WT.cca <- IntegrateLayers(object = WT, method=CCAIntegration,
                                 orig.reduction="pca",new.reduction="integrated.cca",
                                 verbose = FALSE)

integrated_TK.rpca <- IntegrateLayers(object = TK, method=RPCAIntegration,
                                 orig.reduction="pca",new.reduction="integrated.rpca",
                                 verbose = FALSE)

integrated_WT.rpca <- IntegrateLayers(object = WT, method=RPCAIntegration,
                                 orig.reduction="pca",new.reduction="integrated.rpca",
                                 verbose = FALSE)

# 5. find neighbors and clusters and embed using the new integrated assay, for each 'SampleType' individually
# Look at elbow plot
ElbowPlot(integrated_TK.cca)
ElbowPlot(integrated_WT.cca)
ElbowPlot(integrated_TK.rpca)
ElbowPlot(integrated_WT.rpca)


# 6. Do Plots here

Generate_Plots <- function(data,method,dim_num){
  
  data_name <- deparse(substitute(data))
  print(data_name)
  if(str_detect(data_name,'TK')){
    data_type <- 'TK'
  } else{
    data_type <- 'WT'
  }
  if(method == 'RPCA'){
    cluster.name <- 'rpca_clusters'
    reduction.method <- 'integrated.rpca'
    reduction.umap <- 'umap.rpca'
  } else{
    cluster.name <- 'cca_clusters'
    reduction.method <- 'integrated.cca'
    reduction.umap <- 'umap.cca'
  }
  
  # integrated_TK <- FindNeighbors(integrated_TK,reduction='integrated.rpca',dims=1:30)
  # integrated_TK <- FindClusters(integrated_TK,resolution=0.2,cluster.name='rpca_clusters')
  # integrated_TK <- RunUMAP(integrated_TK,reduction='integrated.rpca',dims=1:30,reduction.name='umap.rpca')
  
  
  data <- FindNeighbors(data,reduction=reduction.method,dims=1:dim_num)
  data <- FindClusters(data,resolution=0.2,cluster.name=cluster.name)
  data <- RunUMAP(data,reduction=reduction.method,dims=1:dim_num,reduction.name=reduction.umap)
  
  filename_1 <- paste(data_type,'_Cluster_PC',as.character(dim_num),'_', method,'.png',sep='')
  print(filename_1)
  png(filename_1,width=30,height=10,units='in',res=300)
  p <- DimPlot(data,reduction=reduction.umap,split.by='Condition',group.by=cluster.name,pt.size=0.5,shuffle=T,combine=FALSE)
  print(p)
  dev.off()
  
  filename_2 <- paste('QC_plots_nFeature_RNA_',data_type,'_',method,'_','PC',as.character(dim_num),'.png',sep='')
  print(filename_2)
  png(filename_2,width=30,height=30,units='in',res=300)
  p1 <- VlnPlot(data,c('nFeature_RNA'),pt.size = 0.01,split.by = 'Condition',group.by = cluster.name)
  p2 <- VlnPlot(data,c('nCount_RNA'),pt.size = 0.01,split.by = 'Condition',group.by = cluster.name)
  p3 <- VlnPlot(data,c('percent.mt'),pt.size = 0.01,split.by = 'Condition',group.by = cluster.name)
  print(plot_grid(p1,p2,p3,nrow=3,ncol=1))
  dev.off()
  
  if(method == 'RPCA'){
    dist<- table(data$rpca_clusters)
    write.table(dist,file = paste('dist_',data_type,'_RPCA_',as.character(dim_num),'.csv',sep=''),sep = ',',row.names = T,col.names = NA)
  } else{
    dist<- table(data$cca_clusters)
    write.table(dist,file = paste('dist_',data_type,'_CCA_',as.character(dim_num),'.csv',sep=''),sep = ',',row.names = T,col.names = NA)
  }
  return(data)
}

TK_cca_20 <- Generate_Plots(integrated_TK.cca,'CCA',20)
WT_cca_20 <- Generate_Plots(integrated_WT.cca,'CCA',20)
TK_cca_30 <- Generate_Plots(integrated_TK.cca,'CCA',30)
WT_cca_30 <- Generate_Plots(integrated_WT.cca,'CCA',30)
TK_cca_40 <- Generate_Plots(integrated_TK.cca,'CCA',40)
WT_cca_40 <- Generate_Plots(integrated_WT.cca,'CCA',40)

TK_rpca_20 <- Generate_Plots(integrated_TK.rpca,'RPCA',20)
WT_rpca_20 <- Generate_Plots(integrated_WT.rpca,'RPCA',20)
TK_rpca_30 <- Generate_Plots(integrated_TK.rpca,'RPCA',30)
WT_rpca_30 <- Generate_Plots(integrated_WT.rpca,'RPCA',30)
TK_rpca_40 <- Generate_Plots(integrated_TK.rpca,'RPCA',40)
WT_rpca_40 <- Generate_Plots(integrated_WT.rpca,'RPCA',40)


### Sam made this function to plot distributions......
MetadataDist <- function(object,
                         metadata.1,
                         metadata.2,
                         chunks = 10){
  
  ### TASK 1: Break into bootstrapped samples (without replacement, so no measurement is used more than once)
  # First, let's break the object into the number of chunks
  num.cells.per.chunk <- ceiling(ncol(object)/chunks)
  message(paste('Breaking input object into',chunks,'chunks of approximately',num.cells.per.chunk,'cells each...'))
  # Create a new metadata slot that divvies the input into bootstrapped samples
  bootstrap <- c()
  for(i in 1:chunks){
    bootstrap <- c(bootstrap, rep(paste('Sample',i,sep = "_"),num.cells.per.chunk))
  }
  # remove remaining one cell at the end
  bootstrap <- bootstrap[1:length(bootstrap)-1]
  # check that the length is right
  length(bootstrap) == ncol(object)
  #table(bootstrap)
  # randomize the order completely
  bootstrap <- sample(bootstrap)
  # add to object to use for subsetting
  object$bootstrap <- bootstrap
  split <- SplitObject(object,split.by = 'bootstrap')
  
  ### TASK 2-3: perform requested calculations and store nicely
  dist.data <- data.frame()
  for(i in 1:length(split)){
    
    # create metadata distribution table
    temp <- table(split[[i]]@meta.data[[metadata.1]],split[[i]]@meta.data[[metadata.2]])
    
    # normalize so that all values of metadata.2 for a given metadata.1 add up to 1
    for(j in 1:nrow(temp)){
      temp[j,] <- temp[j,]/rowSums(temp)[j]
    }
    # rowSums(temp) # should all equal 1
    # convert to data frame
    temp <- as.data.frame(temp)
    # add sample ID
    temp$bootstrap <- names(split[i])
    # concatenate
    dist.data <- rbind(dist.data,temp)  
  }
  
  #View(dist.data)
  
  # Provide function outputs
  return(dist.data)
}

### Now let's test it by making some plots. First, let's do Condition vs. clusters(cca/rpca):

Generate_Dist <- function(data,dim_1,dim_2,method,dim_num,data_type){
  dist.data <- MetadataDist(object = data,
                            metadata.1=dim_1,
                            metadata.2=dim_2,
                            chunks = 10)
  
  title_name_1 <- paste(dim_1,' Distribution over ',dim_2,sep='')
  xlab_name_1 <- paste(method,'_cluster',sep='')
  ylab_name_1 <- paste('Percentage of ',dim_1,sep='')
  
  ### make pretty ggplots
  # see https://rpkgs.datanovia.com/ggpubr/reference/geom_pwc.html for stats code
  require(ggplot2)
  require(ggpubr)
  colors.use <- c('#93B7BE','#F19A3E','#3D3B8E','#E072A4',"#B22222","grey","#A034F0",'yellow','#8B786D','green','violet','red','blue','#00BFB2')
  output.plot.1 <-  ggplot(data = dist.data,
                           aes(x = Var2,y=Freq*100,fill=Var1,color=Var1))+
    geom_violin()+
    geom_point(position = position_jitterdodge(dodge.width = 0.9,jitter.width=0.25),size=0.1,color='black')+
    theme_classic()+
    ggtitle(title_name_1)+
    ylab(ylab_name_1)+
    xlab(xlab_name_1)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(face = "bold",hjust = 0.5))+
    geom_pwc(
      aes(group = Var1), 
      remove.bracket = F,
      tip.length = 0,
      method = "t_test",
      label = "{p.adj.signif}",
      y.position = c(37,24,10),
      #label = "{p.adj.format}{p.adj.signif}",
      p.adjust.method = "bonferroni", 
      p.adjust.by = "panel",
      hide.ns = TRUE)+
    guides(fill=guide_legend(title="Condition"))+
    guides(color=guide_legend(title="Condition"))+
    scale_fill_manual(values = colors.use)+
    scale_color_manual(values = colors.use)

  
  ## And if we reverse the ordering, and remake the plot to look good:
  
  dist.data <- MetadataDist(object = data,
                            metadata.1 = dim_2,
                            metadata.2 = dim_1,
                            chunks = 10)
  
  title_name_2 <- paste(dim_2,' Distribution over ',dim_1,sep='')
  xlab_name_2 <- paste(method,'_cluster',sep='')
  ylab_name_2 <- paste('Percentage of ',dim_2)
  
  ### make pretty ggplots
  # see https://rpkgs.datanovia.com/ggpubr/reference/geom_pwc.html for stats code
  require(ggplot2)
  require(ggpubr)
  colors.use <- c('#93B7BE','#F19A3E','#3D3B8E','#E072A4',"#B22222","grey","#A034F0",'yellow','#8B786D','green','violet','red','blue','#00BFB2')
  output.plot.2 <-  ggplot(data = dist.data,
                           aes(x = Var1,y=Freq*100,fill=Var2,color=Var2))+
    geom_violin()+
    geom_point(position = position_jitterdodge(dodge.width = 0.9,jitter.width=0.25),size=0.1,color='black')+
    theme_classic()+
    ggtitle(title_name_2)+
    ylab(ylab_name_2)+
    xlab(xlab_name_2)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(face = "bold",hjust = 0.5))+
    geom_pwc(
      aes(group = Var2), 
      remove.bracket = F,
      tip.length = 0,
      method = "t_test",
      label = "{p.adj.signif}",
      y.position = 80,
      #label = "{p.adj.format}{p.adj.signif}",
      p.adjust.method = "bonferroni", 
      p.adjust.by = "panel",
      hide.ns = T)+
    guides(fill=guide_legend(title="Condition"))+
    guides(color=guide_legend(title="Condition"))+
    scale_fill_manual(values = colors.use)+
    scale_color_manual(values = colors.use)
  
  png(paste(data_type,'_',method,'_PC',as.character(dim_num),'_dist_Plot.png',sep=''),width=35,height=30,units='in',res=300)
  print(plot_grid(output.plot.1,output.plot.2,nrow=2,ncol=1))
  dev.off()
}

Generate_Dist(TK_cca_20,'Condition','cca_clusters','CCA',20,'TK')
Generate_Dist(WT_cca_20,'Condition','cca_clusters','CCA',20,'WT')

Generate_Dist(TK_cca_30,'Condition','cca_clusters','CCA',30,'TK')
Generate_Dist(WT_cca_30,'Condition','cca_clusters','CCA',30,'WT')

Generate_Dist(TK_cca_40,'Condition','cca_clusters','CCA',40,'TK')
Generate_Dist(WT_cca_40,'Condition','cca_clusters','CCA',40,'WT')


Generate_Dist(TK_rpca_20,'Condition','rpca_clusters','RPCA',20,'TK')
Generate_Dist(WT_rpca_20,'Condition','rpca_clusters','RPCA',20,'WT')

Generate_Dist(TK_rpca_30,'Condition','rpca_clusters','RPCA',30,'TK')
Generate_Dist(WT_rpca_30,'Condition','rpca_clusters','RPCA',30,'WT')

Generate_Dist(TK_rpca_40,'Condition','rpca_clusters','RPCA',40,'TK')
Generate_Dist(WT_rpca_40,'Condition','rpca_clusters','RPCA',40,'WT')


### After viewing all the plots, we decided to take method CCA, PC number 40.

### save the data for further use first!

save(TK_cca_40,file = 'TK_cca_40.Robj')
save(WT_cca_40,file = 'WT_cca_40.Robj')

load("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/TK_cca_40.Robj")
load("/Users/nuoyawang/Dropbox/Lung Biology/Yaqing/WT_cca_40.Robj")

### Feature plots for epithelial, endothelium, immune and mesenchyme
png("Feature_Plot_TK_cca_40.png",width=10,height=8,units='in',res=300)
FeaturePlot(TK_cca_40,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)+
  plot_annotation(title ='FeaturePlot TK CCA PC40')&
  theme(plot.title = element_text(size = 20,face = 'bold',hjust = 0.5))
dev.off()

png("Feature_Plot_WT_CCA_40.png",width=10,height=8,units='in',res=300)
FeaturePlot(WT_cca_40,c('Epcam','Cdh5','Col1a1','Ptprc'),label=T)+
  plot_annotation(title ='FeaturePlot WT CCA PC40')&
  theme(plot.title = element_text(size = 20,face = 'bold',hjust = 0.5))
dev.off()

### Generate marker lists
TK_cca_40 <- JoinLayers(TK_cca_40)
mark_TK <- FindAllMarkers(TK_cca_40,only.pos=T)
mark_TK$ratio <- mark_TK$pct.1/mark_TK$pct.2
mark_TK$power <- mark_TK$ratio*mark_TK$avg_log2FC
View(mark_TK)
write.table(mark_TK,file = 'TK_integrated_Markers.csv',sep = ',',row.names = T,col.names = NA)

WT_cca_40 <- JoinLayers(WT_cca_40)
mark_WT <- FindAllMarkers(WT_cca_40,only.pos=T)
mark_WT$ratio <- mark_WT$pct.1/mark_WT$pct.2
mark_WT$power <- mark_WT$ratio*mark_WT$avg_log2FC
View(mark_WT)
write.table(mark_WT,file = 'WT_integrated_Markers.csv',sep = ',',row.names = T,col.names = NA)


