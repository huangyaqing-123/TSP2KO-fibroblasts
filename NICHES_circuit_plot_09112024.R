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
require(patchwork)
require(SeuratWrappers)
require(reticulate)
require(stringr)
require(ggpol)
require(poolr)
require(ggrepel)
require(igraph)

# load data
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/fibro.integrated.Robj')
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/signaling.markers.Robj')
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/Markers_Differentiate_TK.Robj')
load('~/project/scRNA seq analysis 04012024/NICHE_07242024/Markers_Differentiate_WT.Robj')
load('~/project/scRNA seq analysis 04012024/All_Annotated_V2.Robj')

#------ Pdgfb signaling-----
to_plot_1 <- Data_for_CircuitPlot_1(c('Pdgfb—Pdgfrb')) # using the first function for single pair
# Additionally,I exclude Endo, Epi, and Immune in the sending type due to their high expression of PDGFb

# generate weight by normalizing mean globally
to_plot_1$weight <- as.numeric(to_plot_1$mean)/sum(as.numeric(to_plot_1$mean))*100

# split based on genotype
tmp<- split(to_plot_1,to_plot_1$SampleType)
TK<-tmp$TK
WT<-tmp$WT

sending_TK <- TK$SendingType
receiving_TK <- TK$ReceivingType
weight=TK$weight
Pdgfb_TK <- data.frame(from=sending_TK,
                       to=receiving_TK,
                       weight=weight)

sending_WT <- WT$SendingType
receiving_WT <- WT$ReceivingType
weight=WT$weight
Pdgfb_WT <- data.frame(from=sending_WT,
                       to=receiving_WT,
                       weight=weight)

# Load (DIRECTED) graph from data frame 
g <- graph.data.frame(Pdgfb_WT, directed=TRUE)


#----- prepare igraph parameters-----

#add color and labels for the vertex
# define colors
colors.use <- c('#59C9A5','#F18F01',
                '#E072A4','#5FA8D3',
                '#0267C1', '#B3001B','#A034F0','#C490D1', '#7C6A0A','#053B06',
                'green','#53599A','#F0CF65')
names(colors.use) <- unique(fibro.integrated$ReceivingType)

# define labels
# try to apply labels using the below function but did not work so I manually added the label 
#require(expss)
#label <- apply_labels(V(g), 'Dcn+_Fibroblasts'='Dcn+_Fb',
                    # 'Dkk2+_Fibroblasts' = 'Dkk2+_Fb', 
                    # 'Rspo3+_Frzb+_Fibroblasts' = 'Rspo3+_Fb',
                    # 'Fibroblasts_Cycling' = 'Fb_Cycling',
                    #  'Schwann_Cycling' = 'Sc_Cycling',
                    # 'Schwann_Cells' = 'Sc_Precursors',
                    #  'IFN_Stimulated_Fibroblasts' = 'IFN_Stimulated_Fb',
                    # 'Adipogenic_Precursor_Cells' = 'Ebf2+_Fb',
                    #  'Taco1+_Fibroblasts' = 'Taco1+_Fb',
                    #  'Myh11+_Fibroblasts' = 'VSMC',
                    #  'Immune'='Immune',
                    #  'Epithelial' = 'Epithelial',
                    #  'Endothelial'='Endothelial')


label_name_2 <- c('Ebf2+_Fb', 'Dcn+_Fb', 'Dkk2+_Fb','Fb_Cycling','IFN_Stimulated_Fb',
                  'VSMC', 'Rspo3+_Fb','Sc_Precursors', 'Sc_Cycling','Taco1+_Fb',
                  'Endothelial','Epithelial','Immune') # the order of pdgfb nodes
                
V(g)$color <- colors.use 
V(g)$label <- label_name_2
#V(g)$size <- data0.5[1,] # initially I want to use cell proportation to be the size of the nodes,
                          # but then figure looks very weird so I decide no
curves <-autocurve.edges(g)
E(g)$curves <- curves

# make the labels outside
require(scales)
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
n <- 13
lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)

# Plot graph
plot(g, layout = layout.circle(g),
     edge.arrow.size = 0.2,
     vertex.label.cex = 1,
     vertex.label.dist = 4,
     vertex.label.family = "arial", 
     #vertex.label.color = colors.use.2,
     edge.width=weight,
     vertex.label.degree=lab.locs)
  
     
##----BMP signaling pathway total-----
FOI <- c('Bmp4—Acvr2a','Bmp4—Bmpr2', 'Bmp4—Bmpr1a', 'Bmp4—Acvr1')

FOI_single <- c('Bmp4—Bmpr2')

to_plot_1a <- Data_for_CircuitPlot_1('Bmp4—Bmpr1a')
to_plot_2 <- Data_for_CircuitPlot_1('Bmp4—Bmpr2')
to_plot_2a <- Data_for_CircuitPlot_1('Bmp4—Acvr2a')
to_plot_r1 <- Data_for_CircuitPlot_1('Bmp4—Acvr1')

to_plot_bmpr <- cbind(to_plot_1a, to_plot_2, to_plot_2a, to_plot_r1)
to_plot_bmpr$avg_mean <- mean(to_plot_bmpr$weight, to_plot_bmpr$weight)

to_plot_1$strength <- as.numeric(to_plot_1$mean)*as.numeric(to_plot_1$count)
to_plot_1$weight <- as.numeric(to_plot_1$strength)/sum(as.numeric(to_plot_1$strength))*100
to_plot_1$weight.2 <- as.numeric(to_plot_1$mean)/sum(as.numeric(to_plot_1$mean))*100

to_plot_1a$weight.2 <- as.numeric(to_plot_1a$mean)/sum(as.numeric(to_plot_1a$mean))*100
to_plot_2$weight.2 <- as.numeric(to_plot_2$mean)/sum(as.numeric(to_plot_2$mean))*100
to_plot_2a$weight.2 <- as.numeric(to_plot_2a$mean)/sum(as.numeric(to_plot_2a$mean))*100
to_plot_r1$weight.2 <- as.numeric(to_plot_r1$mean)/sum(as.numeric(to_plot_r1$mean))*100

to_plot_1a$weight <- as.numeric(to_plot_1a$mean)
to_plot_2$weight <- as.numeric(to_plot_2$mean)
to_plot_2a$weight <- as.numeric(to_plot_2a$mean)
to_plot_r1$weight <- as.numeric(to_plot_r1$mean)

# weight for edges
#to_plot_bmpr2$weight <- as.numeric(to_plot_bmpr2$mean)
#to_plot_bmpr1$weight <- as.numeric(to_plot_bmpr1$mean)/sum(as.numeric(to_plot_bmpr1$mean))*100

# split into two
tmp<- split(to_plot_r1,to_plot_r1$SampleType)
TK<-tmp$TK
WT<-tmp$WT

sending_TK <- TK$SendingType
receiving_TK <- TK$ReceivingType
weight=TK$weight
Bmp_TK<- data.frame(from=sending_TK,
                    to=receiving_TK,
                    weight=weight)
sending_WT <- WT$SendingType
receiving_WT <- WT$ReceivingType
weight=WT$weight
Bmp_WT <- data.frame(from=sending_WT,
                     to=receiving_WT,
                     weight=weight)

# Load (DIRECTED) graph from data frame 
g_wt <- graph_from_data_frame(Bmp_WT, directed=TRUE)
g_tk <- graph_from_data_frame(Bmp_TK, directed=TRUE)

# define colors and labels
label_name_3 <-  c('Ebf2+_Fb', 'Dcn+_Fb', 'Dkk2+_Fb','Fb_Cycling','IFN_Stimulated_Fb',
                   'VSMC', 'Rspo3+_Fb','Sc_Precursors', 'Sc_Cycling','Taco1+_Fb',
                   'Endothelial','Epithelial',
                   'Immune') # the order of bmp4 nodes

V(g_wt)$color <- colors.use
V(g_wt)$label <- label_name_3
#V(g_wt)$size <- data0.5[1,]
curves <-curve_multiple(g_wt)
E(g_wt)$curves <- curves


# to plot
plot(g_wt, layout = layout.circle(g_wt),
     edge.arrow.size = 0.2,
     vertex.label.cex = 1,
     vertex.label.dist = 4,
     vertex.label.family = "arial", 
     vertex.label.color = colors.use,
     edge.width = weight + 0.00000000000001,
     vertex.label.degree=lab.locs)

# feature plot
FeaturePlot(fibro.integrated,'Bmp4—Bmpr1a', split.by = SampleType)



