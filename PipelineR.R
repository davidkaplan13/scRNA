library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)

setwd('/Users/davidkaplan/Desktop/neuroblastoma_data/kildisuite')
gosh <- readRDS("nb_GOSH.rds")
ag <- readRDS("adr_all.rds")

## Functions
QC_metrics <- function(data){
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  boxplot(data$nFeature_RNA,xlab="nFeature_RNA")
  boxplot(data$nCount_RNA,xlabel="nCount_RNA")
  boxplot(data$percent.mt,xlabel="percent.mt")
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


#HVG
highly_variable_data <- FindVariableFeatures(data_kildisuite, selection.method ="vst", nfeatures=2000)
top10 <- head(VariableFeatures(highly_variable_data),10)
vfp <- VariableFeaturePlot(highly_variable_data)
LabelPoints(plot=vfp, points=top10, repel=TRUE)