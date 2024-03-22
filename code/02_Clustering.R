#########################################################################
##
## RstatClub_2: Multiome Analysis with Seurat and Signac #2
##
## Input: Truncated H5, meta-data.cvs and fragments.tvs
## Output: rds Seurat objects with ATAC assay, one with QCs and other w/o QCs
##.        GEX and ATAC plots
## Clustering Determination
##

## Authors. CSC 
## From: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac
##
########################################################################

# load libraries
library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(here)
set.seed(1234)

here::here()


## Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/GEX_ATAC_preprocessing/"))) {
  dir.create(here("processed-data/GEX_ATAC_preprocessing/"))
}
if (!dir.exists(here("plots/GEX_ATAC_preprocessing/"))) {
  dir.create(here("plots/GEX_ATAC_preprocessing/"))
}


########################    Initials.  ########################  

## Sample: Flash-Frozen Human Healthy Brain Tissue (3k)

## Load pre-existing Seurat object

sample_sample <- 'FFB_Healthy'
# 'FFB_Healthy_QC_QCed.rds'
# 'FFB_Healthy_QC_GEX_ATAC.rds'

messageage('Loading sample: ',sample_sample)

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(sample_sample,'_QC_GEX_ATAC.rds'))
SeuratOBJ <- readRDS(rds_name)

message('Seurat with ATAC object loaded!')   

print(SeuratOBJ)



######## Clustering for RNA

## There are two approaches for running the RNA analysis
##    (1) Standard seurat workflow
##    (2) SCTransform 

## SCTransform model normalization is followed by PCA and UMAP dimensionality reduction
##       Function replaces NormalizeData(), ScaleData(), and FindVariableFeatures()


## (1) Standard seurat workflow

DefaultAssay(SeuratOBJ) <- "RNA"

SeuratOBJ <- NormalizeData(SeuratOBJ, 
                           normalization.method = "LogNormalize", 
                           scale.factor = 10000)
# In the LogNormalize method, Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor

tail(SeuratOBJ[["RNA"]]$data, n=3)

SeuratOBJ <- FindVariableFeatures(SeuratOBJ, 
                                  selection.method = "vst", 
                                  nfeatures = 2000) 
# nfeatures define the top variable features to use
# only used when selection.method is set to 'dispersion' or 'vst'

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SeuratOBJ), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SeuratOBJ, raster=FALSE)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


## Scales and centers features in the dataset
## If variables are provided in vars.to.regress, they are individually regressed against each feature

all.genes <- rownames(SeuratOBJ)
SeuratOBJ <- ScaleData(SeuratOBJ, features = all.genes,
                       vars.to.regress = NULL)




## Perform linear dimensional reduction

# We perform PCA on the scaled data

SeuratOBJ <- RunPCA(SeuratOBJ, features = VariableFeatures(object = SeuratOBJ))
# Outputs a list of genes with the most positive and negative loadings, representing modules of genes that exhibit either correlation across single-cells in the dataset

SeuratOBJ[['pca']]
head(Embeddings(SeuratOBJ, reduction = "pca")[, 1:5])
#head(Stdev(SeuratOBJ, reduction = "pca")[1:5])

##  To visualize both cells and features that define the PCA
DimHeatmap(SeuratOBJ, dims = 1, cells = 300, balanced = TRUE)
# cells and features are ordered according to their PCA scores.


VizDimLoadings(SeuratOBJ, dims = 1:2, reduction = "pca")

#DimPlot(SeuratOBJ)  + NoLegend()

## Determine the ‘dimensionality’ of the dataset
ElbowPlot(SeuratOBJ)


## Cluster the cells

SeuratOBJ <- FindNeighbors(SeuratOBJ, dims = 1:10)

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          resolution = 0.8,
                          algorithm = 1,
                          cluster.name ='C.Lovain')
# Returns a Seurat where the idents have been updated with new cluster info
# latest clustering results will be stored in object metadata under 'seurat_clusters'. 
# Note that 'seurat_clusters' will be overwritten everytime FindClusters is run

table(Idents(SeuratOBJ))

SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain")

colnames(SeuratOBJ@meta.data)
head(SeuratOBJ, n=3)

DimPlot(SeuratOBJ, reduction = "umap.lovain")

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          resolution = 2,
                          algorithm = 1,
                          cluster.name ='C.Lovain.r2')

SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain.2")

## run the non-linear dimensional reduction (UMAP/tSNE)

SeuratOBJ@reductions

DimPlot(SeuratOBJ, reduction = "umap.lovain.2")




########### Next run ATAC analysis. ###########

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth

DefaultAssay(SeuratOBJ) <- "ATAC"

SeuratOBJ <- RunTFIDF(SeuratOBJ,
                      method = 1,
                      scale.factor = 10000)

SeuratOBJ <- FindTopFeatures(SeuratOBJ, 
                             min.cutoff = 'q5',
                             verbose = TRUE)
# Set 'q5':  include 95% most common features as the VariableFeatures.
# Set 10:  include features with >10 total counts in the set of VariableFeatures

SeuratOBJ <- RunSVD(SeuratOBJ)

SeuratOBJ <- RunUMAP(SeuratOBJ, 
                     reduction = 'lsi', 
                     dims = 2:50, 
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")



## Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

SeuratOBJ <- FindMultiModalNeighbors(SeuratOBJ, 
                                     reduction.list = list("pca", "lsi"),
                                     dims.list = list(1:50, 2:50))

SeuratOBJ <- RunUMAP(SeuratOBJ, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          graph.name = "wsnn", 
                          algorithm = 3, verbose = FALSE)

SeuratOBJ@reductions

library('ggplot2')

p1 <- DimPlot(SeuratOBJ, reduction = "umap.lovain.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SeuratOBJ, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SeuratOBJ, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pALL <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

table(Idents(SeuratOBJ))

## Base name to save plots
base_name <- levels(SeuratOBJ$`orig.ident`[1])

png_file <- paste0(base_name,'_CLuster_ALL_mt.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(pALL, filename = png_name, height = 4, width = 8)


## perform sub-clustering on a specific clusterto find additional structure

SeuratOBJ <- FindSubCluster(SeuratOBJ, 
    cluster = 0,
    graph.name = 'wsnn', #Name of graph to use for the clustering algorithm
    subcluster.name = "sub.cluster", #name of sub cluster added in the meta.data
    resolution = 0.5,
    algorithm = 1) #1 = original Louvain algorithm

Idents(SeuratOBJ) <- "sub.cluster"

colnames(SeuratOBJ@meta.data)

table(SeuratOBJ@meta.data[['seurat_clusters']])

table(SeuratOBJ@meta.data[['sub.cluster']])


# add annotations
SeuratOBJ <- RenameIdents(SeuratOBJ, '0_0' = 'CD14-A', '0_1' ='CD14-B')
SeuratOBJ$celltype <- Idents(SeuratOBJ)


p1 <- DimPlot(SeuratOBJ, reduction = "umap.lovain.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SeuratOBJ, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SeuratOBJ, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pALL2 <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

table(Idents(SeuratOBJ))

## Base name to save plots
base_name <- levels(SeuratOBJ$`orig.ident`[1])

png_file <- paste0(base_name,'_CLuster_ALL2_mt.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(pALL, filename = png_name, height = 4, width = 8)


## Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_Clusters_mt.rds'))
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC clusters saved!')  





library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

