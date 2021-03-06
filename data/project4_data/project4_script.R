
library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
setwd("/Volumes/BZ/RG Schier/Group/Blockkurs single cell/blockkurs_sc_projects/")

## Load the dataset
counts <- Read10X_h5("data/project4_data/filtered_feature_bc_matrix.h5")
## Initialize the Seurat object with the raw (non-normalized data)
epi50 <- CreateSeuratObject(
  counts = counts$`Gene Expression`
)
epi50
## QC and selecting cells for further analysis
epi50[["percent.mt"]] <- PercentageFeatureSet(epi50, pattern = "^mt-")
VlnPlot(epi50, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(epi50, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(epi50, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

epi50 <- subset(epi50, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
## Normalizing the data
epi50 <- NormalizeData(epi50, normalization.method = "LogNormalize", scale.factor = 10000)
## Identification of highly variable features (feature selection)
epi50 <- FindVariableFeatures(epi50, selection.method = "vst", nfeatures = 2000)
## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(epi50), 10)
## Plot variable features with and without labels
plot1 <- VariableFeaturePlot(epi50)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
## Scaling the data
all.genes <- rownames(epi50)
epi50 <- ScaleData(epi50, features = all.genes)
## Perform linear dimensional reduction
epi50 <- RunPCA(epi50, features = VariableFeatures(object = epi50))
## Determine the ‘dimensionality’ of the dataset
ElbowPlot(epi50)
## Cluster the cells
epi50 <- FindNeighbors(epi50, dims = 1:10)
epi50 <- FindClusters(epi50, resolution = 0.1)
## Run non-linear dimensional reduction (UMAP/tSNE)
epi50 <- RunUMAP(epi50, dims = 1:10)
DimPlot(epi50, reduction = "umap")

FeaturePlot(epi50, features = c("tbxta")) ##  Mesoderm
FeaturePlot(epi50, features = c("slc26a1")) ## YSL
FeaturePlot(epi50, features = c("tp73")) ## Apoptosis-like
FeaturePlot(epi50, features = c("krt8")) ## EVL
FeaturePlot(epi50, features = c("sox3")) ## Ectoderm

##  Find markers for every cluster compared to all remaining cells, report only the positive ones
epi50.markers <- FindAllMarkers(epi50, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

epi50.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(epi50, features = top10$gene) + NoLegend()

## Assigning cell type identity to clusters
new.cluster.ids <- c("Ectoderm", "Mesoderm", "EVL", "YSL", "Apoptosis-like")
names(new.cluster.ids) <- levels(epi50)
epi50 <- RenameIdents(epi50, new.cluster.ids)
DimPlot(epi50, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Compare expression differences between EVL and YSL
de.genes.evl.ysl <- FindMarkers(
  object = epi50,
  ident.1 = "EVL",
  ident.2 = "YSL",
  min.pct = 0.1,
)

## Find the intersection between metabolically related genes and differentially expressed genes
gene.lipid <-fread("data/project4_data/lipid metabolic process.txt",header=F)  
de.genes.evl.ysl.lipid <- de.genes.evl.ysl[rownames(de.genes.evl.ysl)%in%gene.lipid$V2,]

## Order gene based on expression in EVL, The higher the expression level in EVL, the higher the ranking
de.genes.evl.ysl.lipid <- de.genes.evl.ysl.lipid[order(de.genes.evl.ysl.lipid$avg_log2FC,decreasing=T),]

## Use mesoderm cells as reference 
epi50.new <- subset(x = epi50, idents = c("Mesoderm", "EVL", "YSL"))

## Plot expression patterns in different cell types 
DoHeatmap(epi50.new, features = rownames(de.genes.evl.ysl.lipid)) + NoLegend()

