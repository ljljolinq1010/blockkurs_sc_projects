

wt <- read.csv("GSM3187948_Znf_Wt_Counts.csv", row.names = "X")
mut <- read.csv("GSM3187949_Znf_Mut_Counts.csv", row.names = "X")


mat <- as.matrix(wt)
wt_sparse <- Matrix(mat, sparse = T)

mat <- as.matrix(mut)
mut_sparse <- Matrix(mat, sparse = T)


saveRDS(wt_sparse, "GSM3187948_Znf_Wt_Counts_sparse.rds")

saveRDS(mut_sparse, "GSM3187949_Znf_Mut_Counts_sparse.rds")



pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


wt <- CreateSeuratObject(counts = wt_sparse, project = "wt", min.cells = 3, min.features = 200)
mut <- CreateSeuratObject(counts = mut_sparse, project = "wt", min.cells = 3, min.features = 200)

wt <- NormalizeData(wt, normalization.method = "LogNormalize", scale.factor = 10000)
mut <- NormalizeData(mut, normalization.method = "LogNormalize", scale.factor = 10000)

wt <- FindVariableFeatures(wt, selection.method = "vst", nfeatures = 2000)
mut <- FindVariableFeatures(mut, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(wt), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(wt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(wt)
wt <- ScaleData(wt, features = all.genes)

all.genes <- rownames(mut)
mut <- ScaleData(mut, features = all.genes)



wt <- RunPCA(wt, features = VariableFeatures(object = wt))
mut <- RunPCA(mut, features = VariableFeatures(object = mut))

ElbowPlot(wt)
ElbowPlot(mut)


wt <- FindNeighbors(wt, dims = 1:10)
wt <- FindClusters(wt, resolution = 0.5)

mut <- FindNeighbors(mut, dims = 1:10)
mut <- FindClusters(mut, resolution = 0.5)


wt <- RunUMAP(wt, dims = 1:10)
mut <- RunUMAP(mut, dims = 1:10)



DimPlot(wt, reduction = "umap")
DimPlot(mut, reduction = "umap")


