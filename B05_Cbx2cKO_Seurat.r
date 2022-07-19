setwd("C:/00_RawData/00_Array_and_Seq/21-0920_J15_10X_Cbx2cKO/location_specify_here")

library (Seurat)
library (patchwork)
library (cowplot)
library (dplyr)
library (ggplot2)


# Open data
# Make sure all "barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz" files are present in the specified folders
S01.data <- Read10X(data.dir ="../S01/raw")
S02.data <- Read10X(data.dir ="../S02/raw")
S03.data <- Read10X(data.dir ="../S03/raw")
S04.data <- Read10X(data.dir ="../S04/raw")

# Initialize the Seurat objects, each sample separately
S01 <- CreateSeuratObject(counts = S01.data, project = "S01-cKO", min.cells = 3, min.features = 200)
S02 <- CreateSeuratObject(counts = S02.data, project = "S02-WT", min.cells = 3, min.features = 200)
S03 <- CreateSeuratObject(counts = S03.data, project = "S03-cKO", min.cells = 3, min.features = 200)
S04 <- CreateSeuratObject(counts = S04.data, project = "S04-WT", min.cells = 3, min.features = 200)



# ===QC metrics===

# QC metric: obtain mitochondrial RNA percentage
S01[["percent.mt"]] <- PercentageFeatureSet(S01, pattern = "^mt-")
S02[["percent.mt"]] <- PercentageFeatureSet(S02, pattern = "^mt-")
S03[["percent.mt"]] <- PercentageFeatureSet(S03, pattern = "^mt-")
S04[["percent.mt"]] <- PercentageFeatureSet(S04, pattern = "^mt-")

# QC metric: obtain the percent of the counts that the largest gene occupy
apply(S01@assays$RNA@counts, 2,function(x)(100*max(x))/sum(x)) -> S01$Percent.Largest.Gene
apply(S02@assays$RNA@counts, 2,function(x)(100*max(x))/sum(x)) -> S02$Percent.Largest.Gene
apply(S03@assays$RNA@counts, 2,function(x)(100*max(x))/sum(x)) -> S03$Percent.Largest.Gene
apply(S04@assays$RNA@counts, 2,function(x)(100*max(x))/sum(x)) -> S04$Percent.Largest.Gene

# QC metric: Count/feature ratio. 
# Assumed dou/multiplets if RNA counts are abnormally higher than other cells with similar number of features.
S01$CFRatio <- S01$nCount_RNA/S01$nFeature_RNA
S02$CFRatio <- S02$nCount_RNA/S02$nFeature_RNA
S03$CFRatio <- S03$nCount_RNA/S03$nFeature_RNA
S04$CFRatio <- S04$nCount_RNA/S04$nFeature_RNA



# ===Some example QC plots===

VlnPlot(S01, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0.1, y.max = 20000) 
VlnPlot(S01, features = c("percent.mt"), ncol = 3, pt.size = 0.1) 

# Below plot shows outlier (Abnormally nCount_RNA high cells)
FeatureScatter(S01, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Below plot shows low percent.mt & nFeature_RNA > 2000 show distinct population (likely cells)
FeatureScatter(S01, feature1 = "nFeature_RNA", feature2 = "percent.mt")+ylim(c(0, 50))

# Below plot shows nFeature_RNA > 2000 show lower 'percent.largest.gene'.
FeatureScatter(S01, feature1 = "nFeature_RNA", feature2 = "Percent.Largest.Gene")+ylim(c(0,20))

# Below plot to get an idea of CFRatio distribution
FeatureScatter(S01, feature1 = "nFeature_RNA", feature2 = "CFRatio")+ylim(c(0,20))



# ===Multiple filterings steps below===

# Filter with nFeature_RNA bigger than 2000 (see above rationale)
S01 <- subset (S01, subset = nFeature_RNA > 2000)
S02 <- subset (S02, subset = nFeature_RNA > 2000)
S03 <- subset (S03, subset = nFeature_RNA > 2000)
S04 <- subset (S04, subset = nFeature_RNA > 2000)

# Filter with mitochondrial RNA percentage < 7.5%
S01 <- subset (S01, subset = percent.mt < 7.5)
S02 <- subset (S02, subset = percent.mt < 7.5)
S03 <- subset (S03, subset = percent.mt < 7.5)
S04 <- subset (S04, subset = percent.mt < 7.5)

# Filter with Percentage.Largeset.Gene < 3%
S01 <- subset (S01, subset = Percent.Largest.Gene < 3)
S02 <- subset (S02, subset = Percent.Largest.Gene < 3)
S03 <- subset (S03, subset = Percent.Largest.Gene < 3)
S04 <- subset (S04, subset = Percent.Largest.Gene < 3)

# Filter with Count/Feature ratio, rule out outliers
S01 <- subset (S01, subset = CFRatio < 5.5)
S02 <- subset (S02, subset = CFRatio < 5.5)
S03 <- subset (S03, subset = CFRatio < 5.5)
S04 <- subset (S04, subset = CFRatio < 5.5)



# ===Finished filtering cell, start data processing here===

# Global-scaling normalize data (normalization.method = "LogNormalize", scale.factor = 10000)
S01 <- NormalizeData(S01)
S02 <- NormalizeData(S02)
S03 <- NormalizeData(S03)
S04 <- NormalizeData(S04)

# "Calculate a subset of features that exhibit high cell-to-cell variation in the dataset."
S01 <- FindVariableFeatures(S01, selection.method = "vst", nfeatures = 2000)
S02 <- FindVariableFeatures(S02, selection.method = "vst", nfeatures = 2000)
S03 <- FindVariableFeatures(S03, selection.method = "vst", nfeatures = 2000)
S04 <- FindVariableFeatures(S04, selection.method = "vst", nfeatures = 2000)

# Use the genes common for all 4 samples for integration
# S01: 20710, S02: 21319, S03: 21720, S04: 21634, keep_genes: 20021 elements
keep_genes <- Reduce (intersect, list(rownames(S01), rownames(S02), rownames(S03), rownames(S04)))

# List for integration. Note the variable name intAB is for my own reference.
intAB.list <- list(S01, S02, S03, S04)


# ===Perform integration===

# identify anchors
# output: "retained 3432 anchors"
intAB.anchors <- FindIntegrationAnchors(object.list = intAB.list)

# Produce an 'intergrated' assay
intAB.combined <- IntegrateData(anchorset = intAB.anchors, features.to.integrate = keep_genes)

# Future analysis will be done with corrected "integrated" assay. 
# Later will use "RNA" assay. Note original unmodified "RNA" data still retained. 
DefaultAssay (intAB.combined) <- "integrated"

# Run the standard workflow, scaling, PCA, UMAP, etc.
intAB.combined <- ScaleData(intAB.combined)
intAB.combined <- RunPCA (intAB.combined, npcs=20)
intAB.combined <- RunUMAP (intAB.combined, reduction = "pca", dims = 1:12)
intAB.combined <- FindNeighbors(intAB.combined, reduction="pca", dims=1:12)
intAB.combined <- FindClusters(intAB.combined, resolution=0.1)



# ===Visualization plots====
DimPlot(intAB.combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(intAB.combined, reduction="umap", split.by ="orig.ident")

# note FeaturePlot donw with "integrated", but not "RNA" assay.
FeaturePlot(intAB.combined, features=c("Ccnd2"), cols = c("lightgrey", "red"), min.cutoff="q10", max.cutoff="q90")



# ===Cbx2 genotyping analyses===

# Convert default assay to "RNA"
DefaultAssay (intAB.combined) <- "RNA"


# Note, S01-S03_genotypeV3.csv structure
# Barcode,Genotype
# AAAGATGAGCGTTCCG-1_1,WTorHet
# AAAGATGCAGTTTACG-1_1,WTorHet
# .
# AAACGGGAGCCACTAT-1_1,Mut
# .
# AATCCAGAGTGTCCCG-1_1,NotAssigned
# Total 5357 rows including description row.

# Add 'Genotype' metadata to the intAB.combined Seurat object.
S01.S03.genotypeV3 <- read.csv(file = '../01_AmpliconV2/S01-S03_genotypeV3.csv')
rownames (S01.S03.genotypeV3) <- S01.S03.genotypeV3[,1]
intAB.combined <- AddMetaData(object=intAB.combined, metadata = S01.S03.genotypeV3)

# Check UMAPs based on genotypes
DimPlot(intAB.combined, reduction="umap", split.by = "Genotype")



# ===DEG analyses between genotypes===

#Use only population "3", "2", "0", that express Cbx2
intAB.sub320 <- subset (intAB.combined, idents = c("3", "2","0"))
Idents(intAB.sub320) <- intAB.sub320@meta.data$orig.ident

# Set idents to compare between genotypes
Idents(intAB.sub320) <- intAB.sub320@meta.data$Genotype

# Differential expression test (Wilcoxon rank sum, default)
intAB.sub320.markers <- FindMarkers(intAB.sub320, logfc.threshold = 0.1, ident.1="WTorHet", ident.2= "Mut")
intAB.sub320.avg <- as.data.frame (log1p(AverageExpression(intAB.sub320)$RNA))



# ===Below miscellenous commands I used at some point===



# ===Some commands used to make figures===
DimPlot(intAB.combined, reduction="umap", split.by ="orig.ident")
FeaturePlot(intAB.combined, features=c("Ccnd2"), cols = c("lightgrey", "red"), min.cutoff="q10", max.cutoff="q90")
DimPlot(intAB.combined, reduction="umap", split.by = "Genotype")



# ===Summary counts===
counts_per_cell_S01 <- Matrix::colSums(S01)
counts_per_gene_S01 <- Matrix::rowSums(S01)

hist(log10(counts_per_cell_S01+1))



# ===Data import & exporting functions===
Cbx2_intAB_scaled <- FetchData(intAB.combined, var = "Cbx2")
All_intAB <- FetchData(intAB.combined, vars = keep_genes)
write.csv (intAB.sub320.markers, "./intAB.sub320.markers.csv")



# ===Finding number of elements in the clusters===
intAB.count_table <- table(intAB.combined@meta.data$seurat_clusters, intAB.combined@meta.data$orig.ident, intAB.combined@meta.data$Genotype)



# ===Miscellaneous plots

# Histogram of gene expression w/ counts
ggplot(mapping = aes(intAB.combined@assays$RNA@data["Cbx2",])) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
  ggtitle("Expression")

# Dotplots
DotPlot(object=intAB.combined, features = c("Blm","Cbx2"), cols = c("lightgrey", "blue", "lightgrey","blue"), dot.scale = 8, split.by = "orig.ident")+RotatedAxis() 


