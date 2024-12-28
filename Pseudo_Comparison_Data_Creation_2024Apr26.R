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
library(ggpol)
library(poolr)

# Load data
load('/Volumes/T7/Yaqing/All_Annotated_V2.Robj')

# Set Idents
Idents(All_Annotated_V2) <- 'Global_CellType_V2'

Cluster_List <- rownames(table(All_Annotated_V2@meta.data$Global_CellType_V2))


######## Experiment with different distance metrics...... #########

# This is a function calculates distance matrix of two groups
myEuclid <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow=dim(points1)[1], ncol=dim(points2)[1])
  for(i in 1:nrow(points2)) {
    distanceMatrix[,i] <- sqrt(rowSums(t(t(points1)-points2[i,])^2))
  }
  rownames(distanceMatrix) <- rownames(points1)
  colnames(distanceMatrix) <- rownames(points2)
  distanceMatrix
}

myManhattan <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow=dim(points1)[1], ncol=dim(points2)[1])
  for(i in 1:nrow(points2)) {
    distanceMatrix[,i] <- rowSums(abs(t(t(points1)-points2[i,])))
  }
  rownames(distanceMatrix) <- rownames(points1)
  colnames(distanceMatrix) <- rownames(points2)
  distanceMatrix
}

# myCosine <- function(points1,points2){
#   distanceMatrix <- matrix(NA, nrow=dim(points1)[1], ncol=dim(points2)[1])
#   for(i in 1:nrow(points1)){
#     for(j in 1:nrow(points2)) {
#       distanceMatrix[i,j] <- cosine(points1[i,],points2[j,])
#     }
#     
#   }
#   
#   rownames(distanceMatrix) <- rownames(points1)
#   colnames(distanceMatrix) <- rownames(points2)
#   distanceMatrix
#   
# }

# This is the function calculates the centroid coordinate of a given group
myCentroid <- function(points){
  centroid <- matrix(NA,nrow=1,ncol=dim(points)[2])
  for(i in 1:ncol(points)){
    centroid[,i] <- mean(points[,i])
  }
  rownames(centroid) <- 'centroid'
  colnames(centroid) <- colnames(points)
  centroid
}


# This is the function we used to get the differentiation marker list of a dataset
# Notice that if we want to differentiate by other idents, change "Idents(x) <- x$SampleType" to other idents
get_Marker_List <- function(x,min_pct,logfc_threshold){
  DefaultAssay(x) <- 'RNA'
  Idents(x) <- x$SampleType
  mark <- FindAllMarkers(x,min.pct = min_pct,logfc.threshold = logfc_threshold,only.pos = T)
  mark$ratio <- mark$pct.1/mark$pct.2
  mark$power <- mark$ratio*mark$avg_log2FC
  return(mark)
}

# This function is used to implement small clusters so that we have sufficient data points
Implement_Data <- function(all_data,implement_size_factor,Cluster_Name){
  
  obj <- subset(All_Annotated_V2,idents=Cluster_Name)
  cell_num <- dim(obj@meta.data)[1]
  new_cell_num <- cell_num*implement_size_factor

  # Try to get embedding information for the whole data set
  tmp_neighbor <- FindNeighbors(all_data,reduction='integrated.cca',dims=1:40)

  # This is the embedding matrix for all cells
  dist_matrix_all <- Embeddings(obj=tmp_neighbor[['integrated.cca']])[,1:40]

  # This is the embedding matrix for cells in the small cluster
  tmp_neighbor_sub <- subset(tmp_neighbor,idents=Cluster_Name)
  dist_matrix_sub <- Embeddings(obj=tmp_neighbor_sub[['integrated.cca']])[,1:40]



  # dist_with_rank is a Nx1 matrix, rownames=cells, colunmn value=cell distance to the centroid, in the latent space
  sub_centroid <- myCentroid(dist_matrix_sub)
  all_dist <- myManhattan(dist_matrix_all,sub_centroid)
  dist_with_rank <- all_dist[order(all_dist[,1]),,drop=FALSE]
  selected_cells <- rownames(dist_with_rank)[1:new_cell_num]

  index <- colnames(all_data) %in% selected_cells
  new_obj <- all_data[,index]
  table(new_obj$Condition)
  return(new_obj)
}


# This is the function to find self-defined marker list.
# least_sample_num: at least this number of samples contain the marker
# cluster_ID: a character variable defining which cluster we are looking into

Experiment_Marker_List <- function(Pseudo_Marker_List,least_sample_num,cluster_ID){
  
  Sample_Comb <- t(combn(Sample_Name_List,least_sample_num)) # create all possible sample combinations
  comb_num <- length(Sample_Comb)/least_sample_num # get number of combinations
  
  cur_list <- c() # initializing......
  
  # For loop to iterate all combinations to create a list
  for(i in 1:comb_num){
    first_Sample_Name <- Sample_Comb[i,1]
    first_list_flag <- Pseudo_Marker_List[[first_Sample_Name]]['cluster']==cluster_ID&Pseudo_Marker_List[[first_Sample_Name]]['p_val_adj']<0.001
    first_Sample <- Pseudo_Marker_List[[first_Sample_Name]][first_list_flag,]
    select_Sample_Gene <- first_Sample$gene # first list for intersection initializing...
    for(j in 2:least_sample_num){
      cur_Sample_Name <- Sample_Comb[i,j]
      list_flag <- Pseudo_Marker_List[[cur_Sample_Name]]['cluster']==cluster_ID&Pseudo_Marker_List[[cur_Sample_Name]]['p_val_adj']<0.001
      cur_Sample <- Pseudo_Marker_List[[cur_Sample_Name]][list_flag,]
      cur_Sample_Gene <- cur_Sample$gene
      select_Sample_Gene <- intersect(select_Sample_Gene,cur_Sample_Gene)
    }
    cur_list <- union(cur_list,select_Sample_Gene)  # use union to exclude repetition
  }
  return(cur_list)
}

# Test loop

g <- 'IFN_Stimulated_Fibroblasts'

for(g in Cluster_List){
  
  Cluster_Name <- g
  obj <- subset(All_Annotated_V2,idents=Cluster_Name)
  
  # if(min(table(obj$Condition))<10){
  #   obj <- Implement_Data(all_data = All_Annotated_V2,implement_size_factor = 2, Cluster_Name)
  #   print('Implementation Performed')
  # }
  obj_list <- SplitObject(obj, split.by = "SampleType")
  obj_list_TK <- SplitObject(obj_list$Tsp2.Knockout,split.by='Condition')
  obj_list_WT <- SplitObject(obj_list$WildType,split.by='Condition')
  
  Pseudo_Marker_List <- c()
  
  for(TK_pointer in obj_list_TK){
    TK_name <- rownames(table(TK_pointer$Condition))
    for(WT_pointer in obj_list_WT){
      WT_name <- rownames(table(WT_pointer$Condition))
      cur_comb_name <- paste0(TK_name,'_',WT_name,'_pseudo')
      if(length(colnames(TK_pointer))<10){
        TK_pointer <- Implement_Data(all_data=All_Annotated_V2,implement_size_factor = 3,Cluster_Name)
        print(paste0('Performed Implementation for ',Cluster_Name,' ',TK_name))
      }
      if(length(colnames(WT_pointer))<10){
        WT_pointer <- Implement_Data(all_data=All_Annotated_V2,implement_size_factor = 3,Cluster_Name)
        print(paste0('Performed Implementation for ',Cluster_Name,' ',WT_name))
      }
      cur_comb_obj <- merge(TK_pointer,WT_pointer)
      cur_comb_obj <- JoinLayers(cur_comb_obj)
      cur_comb_obj$pseudo_name <- cur_comb_name
      cur_Marker_List <- c()
      cur_Marker_List[[cur_comb_name]] <- get_Marker_List(cur_comb_obj,min_pct=0.25,logfc_threshold = 0.1)
      Pseudo_Marker_List <- c(Pseudo_Marker_List,cur_Marker_List)
    }
  
  }
  
  obj_name <- paste0(Cluster_Name,'_Pseudo_Marker_List.Robj')
  save(Pseudo_Marker_List,file = obj_name)
  
  print(g)
  
  # load_Name <- paste0('/Volumes/T7/Yaqing/',obj_name)
  # load(load_Name)
  
}

# Calculate Marker List which differentiate TK and WT for each celltype

g <- 'Taco1+_Fibroblasts'
for(g in Cluster_List){
  
  Cluster_Name <- g
  obj_name <- paste0(Cluster_Name,'_Pseudo_Marker_List.Robj')
  load_Name <- paste0('/Volumes/T7/Yaqing/',obj_name)
  load(load_Name)
  Sample_Name_List <- names(Pseudo_Marker_List) # the marker list, names of each sample combination for the current cell type
  sample_num <- length(Sample_Name_List) # the number of samples (pseudo)

  # Markers exist in at least two pseudo samples
  TK_list <- Experiment_Marker_List(Pseudo_Marker_List,least_sample_num=2,cluster_ID='Tsp2.Knockout')
  WT_list <- Experiment_Marker_List(Pseudo_Marker_List,least_sample_num=2,cluster_ID='WildType')


  # TK
  cluster_ID <- 'Tsp2.Knockout'
  result <- data.frame()
  
  # Iterate through the sample list to add more rows
  for (j in 1:sample_num){
    cur_Sample <- Sample_Name_List[j]
    cluster_flag <- sapply(Pseudo_Marker_List[[cur_Sample]]['cluster'],function(x) x==cluster_ID)
    selected_cluster <- Pseudo_Marker_List[[cur_Sample]][cluster_flag,]
    selected_flag <- sapply(selected_cluster$gene, function(x) x %in% TK_list)
    if(length(selected_flag)>0){
      selected_genes <- selected_cluster[selected_flag,]
      if(length(row.names(selected_genes))>0){
        selected_genes$SampleID <- cur_Sample
        result <- rbind(result,selected_genes)
      }
    }
  }
  
  # We only want significant markers, with p-value < 0.001
  p_val_flag <- result['p_val_adj']<0.001
  result <- result[p_val_flag,]
  
  # Calculate the score for each gene
  gene_score <- result %>% count(gene)
  names(gene_score)[names(gene_score) == "n"] <- "Score"
  
  # Add score column to our data frame
  all_with_score <- result %>% 
    left_join(gene_score, by='gene')
  
  # Change NA values to numeric 0
  all_with_score <- all_with_score %>% replace(is.na(.), 0)

  Marker_Diff_SampleType_TK <- all_with_score %>%    
    group_by(gene, cluster) %>%    
    summarize(`Avg pct1` = mean(pct.1),
              `Avg pct2` = mean(pct.2),
              `Avg avg_log2FC` = log2(mean(2^avg_log2FC)),
              `Fisher p value` = fisher(p_val)$p,
              `Fisher p value adj` = fisher(p_val_adj)$p,
              `Avg Score` = mean(Score))
  
  Marker_Diff_SampleType_TK$ratio <- Marker_Diff_SampleType_TK$`Avg pct1`/Marker_Diff_SampleType_TK$`Avg pct2`
  Marker_Diff_SampleType_TK$power <- Marker_Diff_SampleType_TK$ratio*Marker_Diff_SampleType_TK$`Avg avg_log2FC`
  
  # WT
  cluster_ID <- 'WildType'
  result <- data.frame()
  
  # Iterate through the sample list to add more rows
  for (j in 1:sample_num){
    cur_Sample <- Sample_Name_List[j]
    cluster_flag <- sapply(Pseudo_Marker_List[[cur_Sample]]['cluster'],function(x) x==cluster_ID)
    selected_cluster <- Pseudo_Marker_List[[cur_Sample]][cluster_flag,]
    selected_flag <- sapply(selected_cluster$gene, function(x) x %in% WT_list)
    if(length(selected_flag)>0){
      selected_genes <- selected_cluster[selected_flag,]
      if(length(row.names(selected_genes))>0){
        selected_genes$SampleID <- cur_Sample
        result <- rbind(result,selected_genes)
      }
    }
  }
  
  # We only want significant markers, with p-value < 0.001
  p_val_flag <- result['p_val_adj']<0.001
  result <- result[p_val_flag,]
  
  # Calculate the score for each gene
  gene_score <- result %>% count(gene)
  names(gene_score)[names(gene_score) == "n"] <- "Score"
  
  # Add score column to our data frame
  all_with_score <- result %>% 
    left_join(gene_score, by='gene')
  
  # Change NA values to numeric 0
  all_with_score <- all_with_score %>% replace(is.na(.), 0)
  
  Marker_Diff_SampleType_WT <- all_with_score %>%    
    group_by(gene, cluster) %>%    
    summarize(`Avg pct1` = mean(pct.1),
              `Avg pct2` = mean(pct.2),
              `Avg avg_log2FC` = log2(mean(2^avg_log2FC)),
              `Fisher p value` = fisher(p_val)$p,
              `Fisher p value adj` = fisher(p_val_adj)$p,
              `Avg Score` = mean(Score))
  
  Marker_Diff_SampleType_WT$ratio <- Marker_Diff_SampleType_WT$`Avg pct1`/Marker_Diff_SampleType_WT$`Avg pct2`
  Marker_Diff_SampleType_WT$power <- Marker_Diff_SampleType_WT$ratio*Marker_Diff_SampleType_WT$`Avg avg_log2FC`
  
  # Assign output names and save the R objects

  obj_TK_name <- paste0('Markers_Differentiate_TK_',Cluster_Name,'.Robj')
  save(Marker_Diff_SampleType_TK,file = obj_TK_name)
  
  obj_WT_name <- paste0('Markers_Differentiate_WT_',Cluster_Name,'.Robj')
  save(Marker_Diff_SampleType_WT,file = obj_WT_name)
  
}


#### Another Test ####

obj_list <- SplitObject(All_Annotated_V2, split.by = "SampleType")
obj_list_TK <- SplitObject(obj_list$Tsp2.Knockout,split.by='Condition')
obj_list_WT <- SplitObject(obj_list$WildType,split.by='Condition')

Pseudo_Marker_List <- c()

for(TK_pointer in obj_list_TK){
  TK_name <- rownames(table(TK_pointer$Condition))
  for(WT_pointer in obj_list_WT){
    WT_name <- rownames(table(WT_pointer$Condition))
    cur_comb_name <- paste0(TK_name,'_',WT_name,'_pseudo')
    cur_comb_obj <- merge(TK_pointer,WT_pointer)
    cur_comb_obj <- JoinLayers(cur_comb_obj)
    cur_comb_obj$pseudo_name <- cur_comb_name
    cur_Marker_List <- c()
    cur_Marker_List[[cur_comb_name]] <- get_Marker_List(cur_comb_obj,min_pct=0.05,logfc_threshold = 0.1)
    Pseudo_Marker_List <- c(Pseudo_Marker_List,cur_Marker_List)
  }
  
}



Sample_Name_List <- names(Pseudo_Marker_List) # the marker list, names of each sample combination for the current cell type
sample_num <- length(Sample_Name_List) # the number of samples (pseudo)

# Markers exist in at least two pseudo samples
TK_list <- Experiment_Marker_List(Pseudo_Marker_List,least_sample_num=2,cluster_ID='Tsp2.Knockout')
WT_list <- Experiment_Marker_List(Pseudo_Marker_List,least_sample_num=2,cluster_ID='WildType')


# Localized
cluster_ID <- 'Tsp2.Knockout'
result <- data.frame()

# Iterate through the sample list to add more rows
for (j in 1:sample_num){
  cur_Sample <- Sample_Name_List[j]
  cluster_flag <- sapply(Pseudo_Marker_List[[cur_Sample]]['cluster'],function(x) x==cluster_ID)
  selected_cluster <- Pseudo_Marker_List[[cur_Sample]][cluster_flag,]
  selected_flag <- sapply(selected_cluster$gene, function(x) x %in% TK_list)
  if(length(selected_flag)>0){
    selected_genes <- selected_cluster[selected_flag,]
    selected_genes$SampleID <- cur_Sample
    result <- rbind(result,selected_genes)
  }
}

# We only want significant markers, with p-value < 0.01
p_val_flag <- result['p_val_adj']<0.001
result <- result[p_val_flag,]

# Calculate the score for each gene
gene_score <- result %>% count(gene)
names(gene_score)[names(gene_score) == "n"] <- "Score"

# Add score column to our data frame
all_with_score <- result %>% 
  left_join(gene_score, by='gene')

# Change NA values to numeric 0
all_with_score <- all_with_score %>% replace(is.na(.), 0)

Marker_Diff_SampleType_TK <- all_with_score %>%    
  group_by(gene, cluster) %>%    
  summarize(`Avg pct1` = mean(pct.1),
            `Avg pct2` = mean(pct.2),
            `Avg avg_log2FC` = log2(mean(2^avg_log2FC)),
            `Fisher p value` = fisher(p_val)$p,
            `Fisher p value adj` = fisher(p_val_adj)$p,
            `Avg Score` = mean(Score))

Marker_Diff_SampleType_TK$ratio <- Marker_Diff_SampleType_TK$`Avg pct1`/Marker_Diff_SampleType_TK$`Avg pct2`
Marker_Diff_SampleType_TK$power <- Marker_Diff_SampleType_TK$ratio*Marker_Diff_SampleType_TK$`Avg avg_log2FC`

save(Marker_Diff_SampleType_TK, file='Marker_Diff_SampleType_TK.Robj')

# Metastatic
cluster_ID <- 'WildType'
result <- data.frame()

# Iterate through the sample list to add more rows
for (j in 1:sample_num){
  cur_Sample <- Sample_Name_List[j]
  cluster_flag <- sapply(Pseudo_Marker_List[[cur_Sample]]['cluster'],function(x) x==cluster_ID)
  selected_cluster <- Pseudo_Marker_List[[cur_Sample]][cluster_flag,]
  selected_flag <- sapply(selected_cluster$gene, function(x) x %in% WT_list) # modify this %in% list!
  if(length(selected_flag)>0){
    selected_genes <- selected_cluster[selected_flag,]
    selected_genes$SampleID <- cur_Sample
    result <- rbind(result,selected_genes)
  }
}

# We only want significant markers, with p-value < 0.01
p_val_flag <- result['p_val_adj']<0.001
result <- result[p_val_flag,]

# Calculate the score for each gene
gene_score <- result %>% count(gene)
names(gene_score)[names(gene_score) == "n"] <- "Score"

# Add score column to our data frame
all_with_score <- result %>% 
  left_join(gene_score, by='gene')

# Change NA values to numeric 0
all_with_score <- all_with_score %>% replace(is.na(.), 0)

Marker_Diff_SampleType_WT <- all_with_score %>%    
  group_by(gene, cluster) %>%    
  summarize(`Avg pct1` = mean(pct.1),
            `Avg pct2` = mean(pct.2),
            `Avg avg_log2FC` = log2(mean(2^avg_log2FC)),
            `Fisher p value` = fisher(p_val)$p,
            `Fisher p value adj` = fisher(p_val_adj)$p,
            `Avg Score` = mean(Score))

Marker_Diff_SampleType_WT$ratio <- Marker_Diff_SampleType_WT$`Avg pct1`/Marker_Diff_SampleType_WT$`Avg pct2`
Marker_Diff_SampleType_WT$power <- Marker_Diff_SampleType_WT$ratio*Marker_Diff_SampleType_WT$`Avg avg_log2FC`

save(Marker_Diff_SampleType_WT, file='Marker_Diff_SampleType_WT.Robj')


