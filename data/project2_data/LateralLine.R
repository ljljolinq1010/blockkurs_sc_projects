library(dplyr)
library(Seurat)
library(patchwork)
library (tidyverse)


# Load and name the dataset
install.packages("hdf5r")
library(hdf5r)

LateralLine.data <- Read10X_h5("Lateral_line_5dpf.h5", use.names = TRUE, unique.features = TRUE)

# Initialize the Seurat object with the raw (non-normalized data).
LateralLine <- CreateSeuratObject(counts = LateralLine.data, project = "LateralLine", min.cells = 3, min.features = 200)

#How many cells in the dataset?
LateralLine

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#Mitochondrial genes start with mt lowercase in zebrafish!
LateralLine[["percent.mt"]] <- PercentageFeatureSet(LateralLine, pattern = "^mt-")

# Visualize QC metrics as a violin plot
#The first violin plot will tell you where the lower bound (background noise) is, the second will hint towards doublets (e.g cells with an unusual number of RNA molecules), and the last one will tell you the quality of the data. 
#For mitochondrial genes, we set up a cutoff of less than 20 percent
VlnPlot(LateralLine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Subsetting after QC with >200 - < 3000 counts and <20% mtRNA
LateralLine <- subset(LateralLine, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
LateralLine

##LHow many cells did you lose?
LateralLine

#Data normalization
LateralLine <- NormalizeData(LateralLine, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features
LateralLine <- FindVariableFeatures(LateralLine, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LateralLine), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(LateralLine)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(LateralLine)
LateralLine <- ScaleData(LateralLine, features = all.genes)

#Dimensionality reduction
LateralLine <- RunPCA(LateralLine, features = VariableFeatures(object = LateralLine))

# Examine and visualize PCA results a few different ways
print(LateralLine[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(LateralLine, dims = 1:2, reduction = "pca")

DimPlot(LateralLine, reduction = "pca")
DimHeatmap(LateralLine, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(LateralLine, dims = 1:15, cells = 500, balanced = TRUE)

#Calculating the dimensionality
#The number of PCA depends on the p-values you get (how informative each one of them is). 
LateralLine <- JackStraw(LateralLine, num.replicate = 100)
LateralLine <- ScoreJackStraw(LateralLine, dims = 1:20)
JackStrawPlot(LateralLine, dims = 1:20)

#Clustering the cells
#Once you finish the pipeline, you can try different number of PCAs and see how this affects clustering.
#You can also check what the resolution is, from 0.2 to 8. This will split the data into clusters, so play a bit with it until you find something that makes sense.
#We will start with 9PCAs and resolution of 1
LateralLine <- FindNeighbors(LateralLine, dims = 1:9)
LateralLine <- FindClusters(LateralLine, resolution = 1)
head(Idents(LateralLine), 5)

#Creating the visualization with UMAP
LateralLine <- RunUMAP(LateralLine, dims = 1:9)
DimPlot(LateralLine, reduction = "umap")

#Looking for genes in the Feature Plot
FeaturePlot(object = LateralLine, features = c("slc2a5"))

#Saving the RDS object
saveRDS(LateralLine, file = "./LateralLine.rds")

#######Extracting Cluster Markers. This will calculate differential gene expression (the genes that explain the most about each cluster), that you can later analyze.

# find markers for every cluster compared to all remaining cells, report only the positive ones
LateralLine.markers <- FindAllMarkers(LateralLine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
LateralLine.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#saving the cluster markers
# Saving file as csv file
write.csv(LateralLine.markers, "./LateralLine_cluster_markers.csv")



