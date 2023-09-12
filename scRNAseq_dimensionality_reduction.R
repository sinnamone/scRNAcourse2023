# ### Load required libraries

library(Seurat)
library(harmony)
library(Matrix)
library(cowplot)
library(ggplot2)

# ### Define workspace 
#

workDir <- './'
dataDir <- 'data'
# dir.create(file.path(workDir, dataDir), recursive = T, showWarnings = FALSE)
file.exists(file.path(workDir, dataDir))

# ### Load data

load(paste0(file.path(workDir, dataDir), '/', 'pbmc_stim.RData'))

ls()

#Create only one Seurat object with all cell
pbmc <- CreateSeuratObject(counts = cbind(stim.sparse, ctrl.sparse), project = "PBMC", min.cells = 5)

#explore seurat object 
# both counts and data contain the raw UMI counts
pbmc 

#sparse count matrix
head(pbmc@assays$RNA@counts) 

#define datasets with the variable stim
pbmc@meta.data$stim <- c(rep("STIM", ncol(stim.sparse)), rep("CTRL", ncol(ctrl.sparse)))

head(pbmc@meta.data)

# ### Normalization

#Log normalize values (using a scaling factor of 10000)
pbmc <- Seurat::NormalizeData(pbmc, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000) #or pbmc <- Seurat::NormalizeData(verbose = FALSE)

#new layer data with log normalized data 
print(pbmc)

#sparse normalized counts matrix stored in data layer 
head(pbmc[["RNA"]]@data) 

# ### Identification of highly variable features (feature selection)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc

pbmc@assays

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

top50 <- head(VariableFeatures(pbmc), 50)
top50

# plot variable features with and without labels
options(repr.plot.height = 6, repr.plot.width = 8)
pdf("./VariableFeaturePlot.pdf")
plot1 <- VariableFeaturePlot(pbmc); plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE); plot2
dev.off()
# ### Scaling the data

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1.
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#
#

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

head(pbmc[["RNA"]]@scale.data, 5) #centered and scaled normalized umi counts are stored in "scale.data"

# ### Perform linear dimensional reduction

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))  #npcs computed = 50

pbmc@reductions

options(repr.plot.height = 5, repr.plot.width = 5)
DimPlot(object = pbmc, reduction = "pca", pt.size = .1)

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

options(repr.plot.height = 7, repr.plot.width = 4*2)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# ### Visualize batch effect 

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "stim", pt.size = .1)
plot_grid(p1,p2)

# ### Batch effect correction

# #### Run Harmony

options(repr.plot.height = 3, repr.plot.width = 6)
pbmc.corrected <- pbmc %>% 
    RunHarmony("stim", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(pbmc.corrected, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc.corrected, reduction = "harmony", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc.corrected, features = "harmony_1", group.by = "stim",  pt.size = .1)
plot_grid(p1,p2)

# +
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "stim") + 
     ggtitle('Not corrected')
p2 <- DimPlot(object = pbmc.corrected, reduction = "harmony", pt.size = .1, group.by = "stim") + 
        ggtitle('Batch correction')

plot_grid(p1,p2)
# -

# ### Determine the ‘dimensionality’ of the dataset using the elbow plot

options(repr.plot.height = 6, repr.plot.width = 6)
ElbowPlot(pbmc,  ndims = 50, reduction = "pca")

# ### Run non-linear dimensional reduction (UMAP/tSNE)

# #### t-SNE visualization not corrected data 

pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:20, perplexity = 30, max_iter = 1000,
    theta = 0.5, eta = 200, num_threads = 0)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(pbmc, reduction = "tsne", group.by = "stim", pt.size=0.8)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(pbmc, reduction = "tsne", group.by = "stim", pt.size=0.8)

# #### UMAP visualization of not corrected data 

pbmc <- RunUMAP(pbmc, dims = 1:20)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(pbmc, reduction = "umap", group.by = "stim", pt.size = .1)

# #### t-SNE corrected data 

pbmc.corrected <- RunTSNE(pbmc.corrected, reduction = "harmony", dims = 1:20, perplexity = 30, max_iter = 1000,
    theta = 0.5, eta = 200, num_threads = 0)

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(pbmc.corrected, reduction = "tsne", group.by = "stim", pt.size=0.8)

# #### UMAP corrected data 

pbmc.corrected <- RunUMAP(pbmc.corrected, reduction = "harmony", dims = 1:20) 

options(repr.plot.height = 8, repr.plot.width = 8)
DimPlot(pbmc.corrected, reduction = "umap", group.by = "stim", pt.size = .1)

# ### UMAP visualization of batch effect before and after correction

# +
options(repr.plot.height = 6, repr.plot.width = 6*2)

p1 <- DimPlot(object = pbmc, reduction = "umap", pt.size = .1, group.by = "stim") + 
     ggtitle('Not corrected')
p2 <- DimPlot(object = pbmc.corrected, reduction = "umap", pt.size = .1, group.by = "stim") + 
        ggtitle('Batch correction')

plot_grid(p1,p2)

# +
options(repr.plot.height = 6, repr.plot.width = 6*2)

p1 <- DimPlot(object = pbmc, reduction = "tsne", pt.size = .1, group.by = "stim") + 
     ggtitle('Not corrected')
p2 <- DimPlot(object = pbmc.corrected, reduction = "tsne", pt.size = .1, group.by = "stim") + 
        ggtitle('Batch correction')

plot_grid(p1,p2)
# -

# ### Visualize features expression 

genes <- list("T cells" = c("CD3D"),
              "CD4 T cells" = c("CD3D", "CD4"), 
              "CD8 T cells" = c("CD3D", "CD8A", "CD8B"), 
              "Naive CD4+ Tcells" = c("SELL", "CCR7"), 
              "NK cells" = c("GNLY", "NKG7"),
              "B cells" = c("CD79A", "MS4A1"),
              "CD14+ Mono" = c("CD14", "LYZ"),
              "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"), 
              "DC" = c("FCER1A", "CST3"),
              "Platelet" = c("PPBP"))

# +
options(repr.plot.height = 5*5, repr.plot.width = 5*4)

FeaturePlot(pbmc.corrected, 
            features = unlist(genes),
            reduction = "umap",
            min.cutoff = "q9")

# +
options(repr.plot.height = 5*3, repr.plot.width = 5*2)

FeaturePlot(pbmc.corrected, 
            features = c("CD3D", "GNLY", "IFI6"), 
            split.by = "stim", 
            max.cutoff = 3, 
            reduction = "umap",
            cols = c("grey", "red"))
# -

sessionInfo()
