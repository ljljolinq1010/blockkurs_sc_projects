library(dplyr)
library(Seurat)
library(patchwork)

# Set working directory and load Seurat object

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k.rds")

Idents(hypo) <- "Subtype"
idents <- levels(hypo@meta.data$Subtype)


# Subset Seurat object to just immune cells, and just Pachon and Surface, and downsample surface

immune <- subset(hypo, idents = c("Erythrocytes","Tcells","Bcells","Mast_cells","Neutrophils","Macrophages","Microglia"))

Idents(immune) <- "morph"
unique(immune@meta.data$morph)

immune.sub <- subset(immune, idents = c("Choy_surface", "Pachon_cave"))

# There were ~2.17232x more cells from Choy_surface than Pachon_cave, so downsample immune cells by same factor

immune.pachon <- subset(immune, idents = c("Pachon_cave"))
immune.choy <- subset(immune, idents = c("Choy_surface"))


immune.choy.ds <- immune.choy[, sample(colnames(immune.choy), size = round(7144/2.17232), replace=F)]

# Extra raw data, and meta data separately for each

immune.pachon.data <- GetAssayData(immune.pachon, slot = "counts")
immune.choy.data <- GetAssayData(immune.choy.ds, slot = "counts")

immune.pachon.meta <- immune.pachon@meta.data[,c("species", "morph", "sex", "Subtype")]
immune.choy.meta <- immune.choy.ds@meta.data[,c("species", "morph", "sex", "Subtype")]

setwd("/Volumes/BZ/RG Schier/Group/Blockkurs single cell/blockkurs_sc_projects/data/project3_data/")

saveRDS(immune.pachon.data, file = "immune_cells_PachonCave.rds")
saveRDS(immune.choy.data, file = "immune_cells_ChoySurface.rds")

write.csv(immune.pachon.meta, file = "immune_PachonCave_meta_data.csv")
write.csv(immune.choy.meta, file = "immune_ChoySurface_meta_data.csv")