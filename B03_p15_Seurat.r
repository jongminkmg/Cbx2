library(dplyr)
library(Seurat)
library(patchwork)

# Make sure all three following files are present in the folder
# matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz

p15.data <- Read10X(data.dir = "./")
p15 <- CreateSeuratObject(counts = p15.data, project = "p15sg", min.cells = 3, min.features = 200)
p15[["percent.mt"]]<-PercentageFeatureSet(p15, pattern="^MT-")

# Obtained a subset with nFeature_RNA bigger than 1000 (as done in Ernst et al., 2019)
# https://pubmed.ncbi.nlm.nih.gov/30890697/
# https://github.com/MarioniLab/Spermatogenesis2018/blob/master/Preprocessing/10X_scRNAseq/Filtering.Rmd
# Note, already done Count (or UMI#) cut by Cellranger
p15sub <- subset(p15, subset = nFeature_RNA > 1000)

# log normalization
p15 <- NormalizeData(p15sub, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
p15 <- FindVariableFeatures(p15, selection.method = "vst", nfeatures = 2000)

# scaled data
all.genes <- rownames(p15)
p15 <- ScaleData(p15, features = all.genes)

# Linear dimensionality reduction (PCA)
p15 <- RunPCA (p15, features = VariableFeatures(object = p15))

# Cluster cells
p15 <- FindNeighbors(p15, dims = 1:15)
p15 <- FindClusters (p15, resolution = 1.5)
p15 <- RunUMAP(p15, dims=1:15)
DimPlot(p15, reduction = "umap", label = TRUE)

# Check cluster identities by marker genes such as Foxc2, Stra8, cKit, Sycp3, etc
FeaturePlot(p15, features=c("Stra8"))

# Obtain subset of spermatogonia and early spermatocytes
# idents number will be different when re-run
p15sg <- subset(p15, idents = c("5", "7", "3", "0", "9","15"))

# Re cluster spermatogonia and early spermatocytes
p15sg <- FindVariableFeatures (p15sg, selection.method = "vst", nfeatures = 2000)
all.genes <-rownames(p15sg)
p15sg <- ScaleData(p15sg, features = all.genes)
p15sg <- RunPCA (p15sg, features = VariableFeatures (object = p15sg))
p15sg <- FindNeighbors(p15sg, dims = 1:10)
p15sg <- FindClusters (p15sg, resolution = 0.25)
p15sg <- RunUMAP(p15sg, dims=1:10)

# Plot to check
DimPlot(p15sg, reduction = "umap", label = TRUE)
FeaturePlot(p15sg, features=c("Stra8"))
