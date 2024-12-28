# Set WD
setwd("~/project/scRNA seq analysis 04012024")

# Set seed
set.seed(2)

# Pakages
require(Seurat)
require(sp)
require(SeuratObject)
require(ggplot2)
require(pheatmap)

#load data 

load("~/project/scRNA seq analysis 04012024/All_Annotated_V2.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Adipogenic_Precursor_Cells.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Dcn+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Dkk2+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Fibroblasts_Cycling.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_IFN_Stimulated_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Myh11+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Rspo3+_Frzb+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Schwann_Cells.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Schwann_Cycling.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Taco1+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Adipogenic_Precursor_Cells.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Dcn+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Dkk2+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Fibroblasts_Cycling.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_IFN_Stimulated_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Myh11+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Rspo3+_Frzb+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Schwann_Cells.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Schwann_Cycling.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Taco1+_Fibroblasts.Robj")


# transfer object file to readable
saveRDS("Marker_Diff_SampleType_TK.Robj", file = "TK marker list", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

readRDS(file, refhook = NULL)
infoRDS(file)

## Global Analysis 
# load data
load("~/project/scRNA seq analysis 04012024/Marker_Diff_SampleType_TK.Robj")
load("~/project/scRNA seq analysis 04012024/Marker_Diff_SampleType_WT.Robj")


#filter genes
Marker_TK_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 6 & Marker_Diff_SampleType_TK$`Avg avg_log2FC` >= 2,]
Marker_WT_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 6 & Marker_Diff_SampleType_WT$`Avg avg_log2FC` >= 1,]
marker_genes = c(Marker_TK_filtered$gene, Marker_WT_filtered$gene)

# GO analysis
# load required packages
require(clusterProfiler)
require(org.Mm.eg.db)
# all genes
allgenes = Marker_Diff_SampleType_TK$gene
allentrez = bitr(allgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[,2]

# TK
TK_genes = Marker_TK_filtered$gene
TK_entrez = bitr(upgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[,2]
TK_goenrich <- enrichGO(gene = TK_entrez, universe = allentrez, OrgDb = "org.Mm.eg.db", pAdjustMethod = "BH",  ont = 'BP', pvalueCutoff = 1, qvalueCutoff = 1,  minGSSize=NULL, maxGSSize=NULL)
TK_goenrichfig = dotplot(TK_goenrich, color = "qvalue", title = 'Upregulated Pathways in KO',  showCategory=10) 
print(TK_goenrichfig)
# WT
WT_genes = Marker_WT_filtered$gene
WT_entrez = bitr(upgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[,2]
allgenes = Marker_Diff_SampleType_WT$gene
allentrez = bitr(allgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[,2]
WT_goenrich <- enrichGO(gene = WT_entrez, universe = allentrez, OrgDb = "org.Mm.eg.db", pAdjustMethod = "BH",  ont = 'BP', pvalueCutoff = 1, qvalueCutoff = 1,  minGSSize=NULL, maxGSSize=NULL)
WT_goenrichfig = dotplot(WT_goenrich, color = "pvalue", title = 'Upregulated Pathways',  showCategory=10)
print(WT_goenrichfig) # VERY WEIRD RESULTS FOR WT



# Prepare data for volcano plot
Marker_TK_filtered <- Marker_TK_filtered %>%
  mutate(log_pval = -log10(`Fisher p value adj`)) %>%
  rownames_to_column(var = "gene")

# Add a column to indicate significant genes (optional)
Marker_TK_filtered <- Marker_TK_filtered %>%
  mutate(significant = ifelse(`Fisher p value adj` < 0.05 & abs('Avg avg_log2FC') > 0.25, "Significant", "Not Significant"))

# Create volcano plot
volcano_plot <- ggplot(marker_genes, aes(x = avg_log2FC, y = log_pval)) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("red", "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme(legend.position = "top")


#save filtered gene list 
write.csv(Marker_TK_filtered, file = "Marker_TK_filtered.csv")

# calculate the mean of genes in the data
Idents(All_Annotated_V2) = All_Annotated_V2@meta.data$Condition
meanmat = AverageExpression(object = All_Annotated_V2)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat = t(as.matrix(meanmat[match(marker_genes,rownames(meanmat)),]))
drawmat = apply(drawmat,2,scale)
rownames(drawmat) = colnames(meanmat)


# Filter out those non-significant genes (- rik)
drawmat =  drawmat[,-grep(pattern = '.*Rik$',colnames(drawmat))]
drawmat =  drawmat[,-grep(pattern = '^Gm.*',colnames(drawmat))]
pheatmap(drawmat, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("#2066a8",  # Dark Blue
                                               "#8ec1da",  # Medium Blue
                                               "#cde1ec",  # Light Blue
                                               "#ededed",  # Grey
                                               "#f6d6c2",  # Light Red
                                               "#d47264",  # Medium Red
                                               "#ae282c"  ))(100), 
         border_color =NA)

## Cell cluster-based gene differential expression pattern
## Cycling
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Schwann_Cycling.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Schwann_Cycling.Robj")

Marker_TK_Sc_Cycling_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 4 & Marker_Diff_SampleType_TK$`power` >= 2.5,]
Marker_WT_Sc_Cycling_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 4 & Marker_Diff_SampleType_WT$`power` >= 0.5,]
marker_Sc_Cycling_genes = c(Marker_TK_Sc_Cycling_filtered$gene, Marker_WT_Sc_Cycling_filtered$gene)
save(marker_Sc_Cycling_genes, file = 'marker_Sc_Cycling_genes.Robj')

load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Fibroblasts_Cycling.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Fibroblasts_Cycling.Robj")

Marker_TK_Fb_Cycling_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 5 & Marker_Diff_SampleType_TK$`Avg avg_log2FC` >= 1,]
Marker_WT_Fb_Cycling_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 5 & Marker_Diff_SampleType_WT$`Avg avg_log2FC` >= 1,]
Marker_Fb_Cycling_filtered = Marker_TK_Fb_Cycling_filtered + Marker_WT_Fb_Cycling_filtered
marker_Fb_Cycling_genes = c(Marker_TK_Fb_Cycling_filtered$gene, Marker_WT_Fb_Cycling_filtered$gene)
save(marker_Fb_Cycling_genes, file = 'marker_Fb_Cycling_genes.Robj')

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

## Schwann cell-specific cluster differential analysis

# Load data
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Schwann_Cells.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Schwann_Cells.Robj")

# Filter genes (score >6)
Marker_TK_Schwann_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 6 & Marker_Diff_SampleType_TK$`Avg avg_log2FC` >= 1,]
Marker_WT_Schwann_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 6 & Marker_Diff_SampleType_WT$`Avg avg_log2FC` >= 1,]
marker_Schwann_genes = c(Marker_TK_Schwann_filtered$gene, Marker_WT_Schwann_filtered$gene)
save(marker_Schwann_genes, file = 'marker_Schwann_genes.Robj')

# Subset Schwann_cells
Schwann_cells <- subset(All_Annotated_V2, idents='Schwann_Cells')
Idents(Schwann_cells) = Schwann_cells@meta.data$Condition
meanmat = AverageExpression(object = Schwann_cells)
meanmat = meanmat$RNA

# Make a matrix for the heatmap
drawmat_sc = t(as.matrix(meanmat[match(marker_Schwann_genes,rownames(meanmat)),]))
drawmat_sc = apply(drawmat_sc,2,scale)
rownames(drawmat_sc) = colnames(meanmat)
drawmat_sc =  drawmat_sc[,!grepl(pattern = '.*Rik$',colnames(drawmat_sc))]
drawmat_sc =  drawmat_sc[,!grepl(pattern = '^Gm.*',colnames(drawmat_sc))]
pheatmap(drawmat_sc, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)

## Dcn+ fibroblast
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Dcn+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Dcn+_Fibroblasts.Robj")

Marker_TK_Dcn_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 7 & Marker_Diff_SampleType_TK$`Avg avg_log2FC` >= 1,]
Marker_WT_Dcn_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 7 & Marker_Diff_SampleType_WT$`Avg avg_log2FC` >= 1,]
marker_Dcn_genes = c(Marker_TK_Dcn_filtered$gene, Marker_WT_Dcn_filtered$gene)

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

dcn <- subset(All_Annotated_V2, idents='Dcn+_Fibroblasts')

Idents(dcn) = dcn@meta.data$Condition
meanmat = AverageExpression(object = dcn)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat_2 = t(as.matrix(meanmat[match(marker_Dcn_genes,rownames(meanmat)),]))
drawmat_2 = apply(drawmat_2,2,scale)
rownames(drawmat_2) = colnames(meanmat)
drawmat_2 =  drawmat_2[,!grepl(pattern = '.*Rik$',colnames(drawmat_2))]
drawmat_2 =  drawmat_2[,!grepl(pattern = '^Gm.*',colnames(drawmat_2))]
pheatmap(drawmat_2, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)

# Rspo3+_Frzb+_ fibroblast
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Rspo3+_Frzb+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Rspo3+_Frzb+_Fibroblasts.Robj")

Marker_TK_Rspo3_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 7 & Marker_Diff_SampleType_TK$`power` >= 1,]
Marker_WT_Rspo3_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 7 & Marker_Diff_SampleType_WT$`power` >= 1,]
marker_Rspo3_genes = c(Marker_TK_Rspo3_filtered$gene, Marker_WT_Rspo3_filtered$gene)

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

r <- subset(All_Annotated_V2, idents='Rspo3+_Frzb+_Fibroblasts')

Idents(r) = r@meta.data$Condition
meanmat = AverageExpression(object = r)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat_r = t(as.matrix(meanmat[match(marker_Rspo3_genes,rownames(meanmat)),]))
drawmat_r = apply(drawmat_r,2,scale)
rownames(drawmat_r) = colnames(meanmat)
drawmat_r =  drawmat_r[,!grepl(pattern = '.*Rik$',colnames(drawmat_r))]
drawmat_r =  drawmat_r[,!grepl(pattern = '^Gm.*',colnames(drawmat_r))]
drawmat_r =  drawmat_r[,!grepl(pattern = '^AA.*',colnames(drawmat_r))]

pheatmap(drawmat_r, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)

# Dkk2+_ fibroblast
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Dkk2+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Dkk2+_Fibroblasts.Robj")

Marker_TK_Dkk2_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 6 & Marker_Diff_SampleType_TK$`power` >= 1,]
Marker_WT_Dkk2_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 6 & Marker_Diff_SampleType_WT$`power` >= 1,]
marker_Dkk2_genes = c(Marker_TK_Dkk2_filtered$gene, Marker_WT_Dkk2_filtered$gene)

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

d <- subset(All_Annotated_V2, idents='Dkk2+_Fibroblasts')

Idents(d) = d@meta.data$Condition
meanmat = AverageExpression(object = d)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat_d = t(as.matrix(meanmat[match(marker_Dkk2_genes,rownames(meanmat)),]))
drawmat_d = apply(drawmat_d,2,scale)
rownames(drawmat_d) = colnames(meanmat)
drawmat_d =  drawmat_d[,!grepl(pattern = '.*Rik$',colnames(drawmat_d))]
drawmat_d =  drawmat_d[,!grepl(pattern = '^Gm.*',colnames(drawmat_d))]
pheatmap(drawmat_d, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)

# Taco1+_ fibroblast
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Taco1+_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Taco1+_Fibroblasts.Robj")

Marker_TK_Taco1_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 4 & Marker_Diff_SampleType_TK$`power` >= 1,]
Marker_WT_Taco1_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 4 & Marker_Diff_SampleType_WT$`power` >= 1,]
marker_Taco1_genes = c(Marker_TK_Taco1_filtered$gene, Marker_WT_Taco1_filtered$gene)

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

t <- subset(All_Annotated_V2, idents='Taco1+_Fibroblasts')

Idents(t) = t@meta.data$Condition
meanmat = AverageExpression(object = t)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat_t = t(as.matrix(meanmat[match(marker_Taco1_genes,rownames(meanmat)),]))
drawmat_t = apply(drawmat_t,2,scale)
rownames(drawmat_t) = colnames(meanmat)
drawmat_t =  drawmat_t[,!grepl(pattern = '.*Rik$',colnames(drawmat_t))]
drawmat_t =  drawmat_t[,!grepl(pattern = '^Gm.*',colnames(drawmat_t))]
pheatmap(drawmat_t, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)

# IFN_Stimulated_ fibroblast
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_IFN_Stimulated_Fibroblasts.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_IFN_Stimulated_Fibroblasts.Robj")

Marker_TK_IFN_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 7 & Marker_Diff_SampleType_TK$`power` >= 10,]
Marker_WT_IFN_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 7 & Marker_Diff_SampleType_WT$`power` >= 10,]
marker_IFN_genes = c(Marker_TK_IFN_filtered$gene, Marker_WT_IFN_filtered$gene)

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

t <- subset(All_Annotated_V2, idents='IFN_Stimulated_Fibroblasts')

Idents(t) = t@meta.data$Condition
meanmat = AverageExpression(object = t)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat_i = t(as.matrix(meanmat[match(marker_IFN_genes,rownames(meanmat)),]))
drawmat_i = apply(drawmat_i,2,scale)
rownames(drawmat_i) = colnames(meanmat)
drawmat_i =  drawmat_i[,!grepl(pattern = '.*Rik$',colnames(drawmat_i))]
drawmat_i =  drawmat_i[,!grepl(pattern = '^Gm.*',colnames(drawmat_i))]
pheatmap(drawmat_i, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)

# Adipogenic_Precursor_Cells
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_TK_Adipogenic_Precursor_Cells.Robj")
load("~/project/scRNA seq analysis 04012024/Markers_Differentiate_WT_Adipogenic_Precursor_Cells.Robj")

Marker_TK_AP_filtered = Marker_Diff_SampleType_TK[Marker_Diff_SampleType_TK$`Avg Score` >= 2 & Marker_Diff_SampleType_TK$`power` >= 6.9,]
Marker_WT_AP_filtered = Marker_Diff_SampleType_WT[Marker_Diff_SampleType_WT$`Avg Score` >= 2 & Marker_Diff_SampleType_WT$`power` >= 2,]
marker_AP_genes = c(Marker_TK_AP_filtered$gene, Marker_WT_AP_filtered$gene)

# calculate the mean of genes in the data
# Identify cell cluster of interest (replace 0 with the cluster ID you are interested in)

a <- subset(All_Annotated_V2, idents ='Adipogenic_Precursor_Cells')

Idents(a) = a@meta.data$Condition
meanmat = AverageExpression(object = a)
meanmat = meanmat$RNA

# make a matrix for the heatmap
drawmat_a = t(as.matrix(meanmat[match(marker_AP_genes,rownames(meanmat)),]))
drawmat_a = apply(drawmat_a,2,scale)
rownames(drawmat_a) = colnames(meanmat)
drawmat_a =  drawmat_a[,!grepl(pattern = '.*Rik$',colnames(drawmat_a))]
drawmat_a =  drawmat_a[,!grepl(pattern = '^Gm.*',colnames(drawmat_a))]
pheatmap(drawmat_a, cluster_cols=F, 
         cluster_rows=F,  
         color = colorRampPalette(colors = c("blue", 'white', "red"))(100), 
         border_color =NA)


