library(dplyr)
library(Seurat)
library(patchwork)

##############################################################################################
# In this project, your goal is to analyse single-data of immune cells from the brains  ######
# of surface and cave-morphs of the Mexican tetra, Astyanax mexicanus. There are lots of #####
# differences, do not worry about finding them all - focus on what interests you! Good luck! #
##############################################################################################

### Load and recreate Seurat Objects and combined object

##### IMPORTANT ######
##### IMPORTANT ######
##### IMPORTANT ######

# FIRTST SET YOUR WORKING DIRECTORY TO WHERE YOU DOWNLOADED THE DATA (change the directory in this script!!)

setwd("/Volumes/BZ/RG Schier/Group/Blockkurs single cell/blockkurs_sc_projects/data/project3_data")

# Read in the count matrices
pachon.data <- readRDS("immune_cells_PachonCave.rds")
surface.data <- readRDS("immune_cells_ChoySurface.rds")

# Generate Seurat objects for each dataset

pachon <- CreateSeuratObject(counts = pachon.data, project = "Pachon_cave", min.cells = 3, min.features = 200)
pachon@meta.data$morph <- "Pachon_cave"
surface <- CreateSeuratObject(counts = surface.data, project = "Surface", min.cells = 3, min.features = 200)
surface@meta.data$morph <- "Choy_surface"

# Load and add meta data which contains cell type annotations

pachon.meta.data <- read.csv("immune_PachonCave_meta_data.csv", head = T, row.names = "X")
surface.meta.data <- read.csv("immune_ChoySurface_meta_data.csv", head = T, row.names = "X")

# In this dataset, we also recorded the animals sex, you can add it to the meta data using these commands

pachon <- AddMetaData(pachon, pachon.meta.data$sex, col.name = "sex")
surface <- AddMetaData(surface, surface.meta.data$sex, col.name = "sex")

# Finally, add the Cell-type meta data
pachon <- AddMetaData(pachon, pachon.meta.data$Subtype, col.name = "cell_type")
surface <- AddMetaData(surface, surface.meta.data$Subtype, col.name = "cell_type")

### Now you can merge the two objects together, which will make analysis simplier ###
###         This step is optional, you can keep them separate or together         ###

immune <- merge(pachon, surface)

## Then continue with the normal Seurat steps to Normalize, FindVariableFeatures, PCA analysis, clustering, and visualization

immune <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)
immune <- FindVariableFeatures(immune, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(immune), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(immune)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(immune)
immune <- ScaleData(immune, features = all.genes)

## Run PC analysis

immune <- RunPCA(immune, features = VariableFeatures(object = immune))

ElbowPlot(immune)

# Decide on how many PCs to use in clustering

n.pc <- 10 # Put your number here
  
immune <- FindNeighbors(immune, dims = 1:n.pc)

# Find clusters at a few resolutions if you want (though cell types are already known)
immune <- FindClusters(immune, resolutio = 0.25)
immune <- FindClusters(immune, resolutio = 0.5)
immune <- FindClusters(immune, resolutio = 1)

# Run UMAP for visualization

immune <- RunUMAP(immune, dims = 1:n.pc)

# Some basic plots to look at the data

by_morph <- DimPlot(immune, reduction = "umap", group.by = "morph")
by_cell_type <- DimPlot(immune, reduction = "umap", group.by = "cell_type", label = T)

by_morph + by_cell_type

# Here is code to compare the number of cells in each cluster between morphs
table(immune@meta.data$morph, immune@meta.data$cell_type)

# Some example marker genes to start
DotPlot(immune, group.by = "cell_type", features = c("cd74a", "apoc1", "adam8a", "cxcr4b", "apoeb"))


## PLACE YOUR CODE HERE ##

