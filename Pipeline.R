library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)

setwd('/project/home22/dk1222/neuroblastoma_analysis')
gosh <- readRDS("nb_GOSH.rds")
ag <- readRDS("adr_all.rds")


## Functions
QC_metrics <- function(data){
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  boxplot(data$nFeature_RNA,xlab="nFeature_RNA")
  boxplot(data$nCount_RNA,xlabel="nCount_RNA")
  boxplot(data$percent.mt,xlabel="percent.mt")
  #
  #p <- ggplot(data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
    #geom_point(size = 0.5) +
    #scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
    #scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
    #theme_bw() +
    #scale_color_viridis(
      #name = 'Percent MT\ntranscripts',
      #limits = c(0,1),
      #labels = scales::percent,
      #guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
    #)
  #p
  
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm')+
    guides(color = NULL)
  
  return(data)
}

mito_filtering <- function(data){
  median_mt <- median(data$percent.mt)
  mad_mt <- mad(data$percent.mt)
  print(sd(data$percent.mt))
  quantile_filter <- quantile(data$percent.mt)[4]
  data <- subset(data, subset = percent.mt < quantile_filter)
  return(data)
  
}

## Call Point
gosh <- QC_metrics(gosh)
ag <- QC_metrics(ag)
control <- mito_filtering(ag)


#Merging
control$batch <- sample(x = 'control', size = ncol(x = control), replace = TRUE)
gosh$batch <- sample(x = 'gosh', size = ncol(x = gosh), replace = TRUE)
data_kildisuite <- merge(x=control, y=gosh, merge.data=TRUE)
data_kildisuite <- NormalizeData(data_kildisuite)


#HVG + Labels for top 20 most variable genes
data_kildisuite <- FindVariableFeatures(data_kildisuite, selection.method ="vst", nfeatures=2000)
top20 <- head(VariableFeatures(data_kildisuite),20)
vfp <- VariableFeaturePlot(data_kildisuite)
LabelPoints(plot=vfp, points=top20, repel=TRUE)

variable_genes <- data_kildisuite@assays$RNA@var.features

# Scaling
data_kildisuite <- ScaleData(data_kildisuite, features=variable_genes)

#PCA and Visualisation
data_kildisuite <- RunPCA(data_kildisuite, features = variable_genes)
ElbowPlot(data_kildisuite,ndims=40)

pct <- data_kildisuite[["pca"]]@stdev / sum(data_kildisuite[["pca"]]@stdev) * 100
cumulative_summary <- cumsum(pct)

co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

pcs <- min(co2)
plot_df <- data.frame(pct = pct, 
                      cumu = cumulative_summary, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumulative_summary, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  theme_bw()


#Batch Correction

Harmony_seurat <- RunHarmony(data_kildisuite, group.by.vars= "batch", reduction.save = "harmony")
Harmony_seurat <- FindNeighbors(Harmony_seurat, reduction = "harmony", k.param = 10, dims=1:17)
Harmony_seurat <- FindClusters(Harmony_seurat, resolution=0.1)
Harmony_seurat <- RunUMAP(Harmony_seurat, reduction="harmony", dims=1:17)

DimPlot(Harmony_seurat, reduction='pca')
DimPlot(Harmony_seurat, reduction='umap',label=TRUE)+ NoLegend()

Harmony_seurat.markers <- FindAllMarkers(Harmony_seurat, only.pos=TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Harmony_seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by=avg_log2FC)


cluster_0 <- Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 0),]
#
cluster_1 <- Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 1),]
#HBG2 / HBG1 / HBB - Erythrocyte
cluster_2 <- Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 2),]
cluster_3 <- Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 3),]
#COL1A1 / COL1A2 / COL3A1 / DCNN / LUM - Fibroblasts

cluster_4 <- Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 4),]
#STMN2 / TUBB2B / TUBA1A /PRPH / GAP43 / NPY / MLLT11 / STMN4 - Sympathoblasts
cluster_5 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 5),])[1:5]
cluster_5
#
cluster_6 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 6),])[1:5]
cluster_6

cluster_7 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 7),])[1:5]
cluster_7

cluster_8 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 8),])[1:5]
cluster_8

cluster_9 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 9),])[1:5]
cluster_9

cluster_10 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 10),])[1:5]
cluster_10

cluster_11 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 11),])[1:5]
cluster_11

cluster_12 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 12),])[1:5]
cluster_12

cluster_13 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 13),])[1:5]
cluster_13

cluster_14 <- rownames(Harmony_seurat.markers[which(Harmony_seurat.markers$cluster == 14),])[1:5]
cluster_14
