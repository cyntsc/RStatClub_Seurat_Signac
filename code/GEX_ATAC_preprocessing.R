########################################################################
## RstatClub script
## Create a Seruat Object
## Measure Quality Controls for GEX and ATAC assays from CellRanger-ARC data
## We are using Seurat & Signac packages
##
## Input: Truncated H5, meta-data.cvs and fragments.tvs
## Output: rds Seurat objects and some plots
##
## NOTES: 
## For a ~10k cells it is recommended ~40G free mem to process the TSS() ATAC score.  
## Without storing the base-resolution matrix of integration counts at each site you can use less memory, but does not allow plotting the accessibility profile at the TSS.
## For slurm env: $srun --pty --mem=40GB --x11 bash
########################################################################
# 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# # For Seurat from CRAN repository
# install.packages('Seurat')     # v4 
# # Seurat v5
# # Current available repository on: https://satijalab.org/seurat/articles/install.html
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# # Seurat Old versions
# remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
# 
# # For Signac
# setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
# install.packages("Signac")
# # Development version
# if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")
# devtools::install_github("stuart-lab/signac", ref = "develop")
# 
# # For Rainer J (2017). EnsDb.Hsapiens.v86: Ensembl based annotation package. R package version 2.99.0
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# install.packages("here")
# install.packages("tidyverse")

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#packageVersion("Seurat")
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(here)

here::here()

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/GEX_ATAC_preprocessing/"))) {
    dir.create(here("processed-data/GEX_ATAC_preprocessing/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/GEX_ATAC_preprocessing/"))) {
    dir.create(here("plots/GEX_ATAC_preprocessing/"))
}

# Some descriptive stats for further analysis
get_descriptive_stats_GEX <- function(SeuratO, sample_name, sample_tissue) {    
    
    tab_stats <- get_basic_stats_GEX(SeuratO, sample_tissue)
    message(paste0('Exporting table with QC quantiles for sample ', sample_name))
    s_file_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(sample_name,'_GEX_MITO_stats.csv'))
    write.csv(tab_stats, file=s_file_name, quote=TRUE, row.names=FALSE)
    return(tab_stats)
}


########################    Initials ########################  
# Recommended 40G of free_mem to 3k-10k cells if TSS mtx is required

# Calculate and save some descriptive stats for further analysis
source(here("code", "Seurat_custom_functions.R"))       # Call functions to read paths

base_path <- paste0('processed-data/cellranger-arc')    # Read H5 file
#s_sample <- 'pbmc3k'
s_sample <- 'FFB_Healthy'
s_tissue <- 'human'
message('Processing sample: ',s_sample, ' from ', s_tissue, ' brain tissue.')

if (s_sample==pbmc3k) {
    # Read the filtered barcode matrix, the meta-data and the atac fragment from cellranger-ARC output
    s_featured_bc_mtx_path <- here(base_path, "human_brain_3k_filtered_feature_bc_matrix.h5" ) 
    # NOTE. Since v2.0.0+ is used the same file name for the meta-data file. 
    s_meta_data_path <- here(base_path, "human_brain_3k_per_barcode_metrics.csv")
    s_atac_path <- here(base_path,"human_brain_3k_atac_fragments.tsv.gz")
    
} else {
    
    s_meta_data_path <- here(base_path, "pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv")
    s_featured_bc_mtx_path <- here(base_path, "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5" )  # Read csv
    s_atac_path <- here(base_path,"pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
    
}
    
# filtered bc mtx
mtx <- Read10X_h5(s_featured_bc_mtx_path)             # Returns a sparse matrix with rows and columns labeled
rna_counts <- mtx$`Gene Expression`
print(head(rna_counts, n = 3))

metadata <- read.csv(file = s_meta_data_path, header = TRUE, row.names = 1)
# subset specific fields in the meta data df
meta_tmp = c('atac_peak_region_fragments','atac_fragments')
meta = metadata[meta_tmp]
print(head(meta, n = 3))

SeuratOBJ <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = s_sample,
    meta.data = meta
)

SeuratOBJ
str(SeuratOBJ)
SeuratOBJ@meta.data
head(SeuratOBJ, n=3)
#SeuratOBJ@active.ident[1]

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(s_sample,'.rds'))
saveRDS(SeuratOBJ, file = rds_name)

## Read pre-existing Seurat
## pre_existing_seurat <- '/Users/ccardin2/Documents/RStatClub2/processed-data/GEX_ATAC_preprocessing/'
#Seurat2 <- LoadSeuratRds(here('processed-data/GEX_ATAC_preprocessing', paste0(s_sample,'.rds')))

# Calculate MITO levels
SeuratOBJ$log10GenesPerUMI <- log10(SeuratOBJ$nFeature_RNA) / log10(SeuratOBJ$nCount_RNA)

# Select the correct Symbol for each genome. By default is assumed MT, belonging to humans.
# For mouse genome we use the Symbol: mt-. 
# https://www.ncbi.nlm.nih.gov/nuccore/34538597 

# Calculate MITO GENES from human
SeuratOBJ$percent.mt <- PercentageFeatureSet(SeuratOBJ, pattern = "^MT-")
SeuratOBJ$MTRatio <- SeuratOBJ$percent.mt / 100 
SeuratOBJ$percent.ribo <- PercentageFeatureSet(SeuratOBJ, pattern = "^RP[LS]")

head(SeuratOBJ, n=3)

message('^MITO and ^RIBO levels processed successfully')



########           Get Visualizations for the GEX           ######## 

# Build and plot:
#       UMI, Genes, MITO and RIBO violin plots
#       Number of cells per sample
#       UMI/transcripts per cell
#       Distribution of genes per cell (histogram)

# Call Vln plot from Seurat suite. (ggplot2 wrapped)    
sfeature <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt')

# Violin plot with UMIs, Genes, ^MT and RIBO levels
VlnPlot(SeuratOBJ, features = sfeature, group.by = 'orig.ident')

# Plot Genes and UMIs by density per cell 
df_genes_per_cell <- as.data.frame(SeuratOBJ[[]])
df_genes_per_cell %>%
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +
    geom_vline(xintercept = 300) +
    ylab("Log10(UMIs)") +
    xlab("Gene-counts") +
    ggtitle("Genes density by cell") 

df_genes_per_cell %>%
    ggplot(aes(x=orig.ident, y=(nFeature_RNA), fill=orig.ident)) +
    geom_boxplot(alpha = 0.7) +
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("Log10(nFeature_RNA)") +
    xlab("") +
    ggtitle("Genes distribution by cell")

# Correlation btw Genes/UMIs 
df_genes_per_cell %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
    #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 200, linetype=2) +
    geom_hline(yintercept = 200, linetype=2) +
    facet_wrap(~orig.ident) +
    ylab("log10(nFeature_RNA)") +
    xlab("log10(UMIs)") +
    ggtitle('UMIs/Genes by MT levels')


#tab_stats <- get_descriptive_stats_GEX(SeuratOBJ, s_sample, s_tissue)


########  Create ATAC assay and attach it to Seurat object     ##### 

# Create gene annotations for hg38 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))
seqlevelsStyle(annotations) <- "UCSC"
length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) 
genome(annotations) <- "hg38"
# names(mcols(annotations))
#[1] "tx_id"        "gene_name"    "gene_id"      "gene_biotype" "type"   


message('Processing ATAC for sample ', levels(SeuratOBJ$`orig.ident`[1]))
base_name <- paste0(levels(SeuratOBJ$`orig.ident`[1]),'_QC')

# Create the chromatin assay with annotations and attach it to the Seurat object 
atac_counts <- mtx$Peaks

# Convert a genomic coordinate string to a GRanges object ands assign it to the ATAC counts mtx
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))      
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
message('ATAC counts annotation attached successfully')

# Create chromatin assay
chrom_assay <- CreateChromatinAssay(counts = atac_counts, 
                                    sep = c(":", "-"), 
                                    fragments = s_atac_path, 
                                    annotation = annotations)
chrom_assay

SeuratOBJ[["ATAC"]] <- chrom_assay
# Annotations of the object are set
Annotation(SeuratOBJ[["ATAC"]]) <- annotations

head(SeuratOBJ, n = 3)

message('Chromatin assay attached to seurat object successfully')



# Calculate Nucleosome Signal

DefaultAssay(SeuratOBJ) <- "ATAC"

# Calculate the strength of the nucleosome signal per cell
SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p1_NS <- FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')
#head(SeuratOBJ, n = 3) 

# Calculate the "Transcription Start Site (TSS)" enrichment score
tryCatch( {
    # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
    # 1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
    # 2. Error in `colnames<-`(`*tmp*`, value = seq_len(length.out = region.width) -  : attempt to set 'colnames' ...
    # This is a vague message that would happen if no fragments are found in the set of TSS regions. 
    # You could double-checking that the correct gene annotations is being used or you have a low ATAC quality. 
    
    SeuratOBJ <- TSSEnrichment(SeuratOBJ, fast = FALSE) 
    # Group by cells with TSS enrichment scores in two groups.
    SeuratOBJ$high.tss <- ifelse(SeuratOBJ$TSS.enrichment > 2, 'High', 'Low')
    #colnames(SeuratOBJ@meta.data)
    p1_TSS <- TSSPlot(SeuratOBJ, group.by = 'high.tss') + NoLegend()
    png_file_TSS <- paste0(base_name,'_TSS.png')
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_TSS)
    ggsave(p1_TSS, filename = png_name, height = 4, width = 4)
        
    }
    , error = function(e) { print('An error occurred. Check the annotation or you probably have ATAC quality loss.') } )
    
    # Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$blacklist_fraction <- FractionCountsInRegion(
        object = SeuratOBJ,
        assay = 'ATAC',
        regions = blacklist_hg38
    )
    # Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$pct_reads_in_peaks <- SeuratOBJ$atac_peak_region_fragments / SeuratOBJ$atac_fragments * 100
    SeuratOBJ$blacklist_ratio <- SeuratOBJ$blacklist_fraction / SeuratOBJ$atac_peak_region_fragments
    #        Error in `x[[i, drop = TRUE]]`:
    #       ! 'blacklist_fraction' not found in this Seurat object
    # Plot Peaks in black ratio and ATAC main feature scoreds
    p1_BlackR <-  VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2)  
    p1_ATAC <- VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), group.by = "orig.ident") 
    
    png_file_NS <- paste0(base_name,'_Fragment_Distribution_grp.png')
    png_file_BlackR <- paste0(base_name,'_reads_in_peaks.png')
    png_file_ATAC <- paste0(base_name,'_ATAC_QCs.png')
    
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_NS)
    ggsave(p1_NS, filename = png_name, height = 4, width = 4)
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_BlackR)
    ggsave(p1_BlackR, filename = png_name, height = 4, width = 5)
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_ATAC)
    ggsave(p1_ATAC, filename = png_name, height = 4, width = 7)
    
    message('ATAC QCs plots saved!')  
    

# Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_GEX_ATAC.rds'))
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC object saved!')   



########  Standard thresholds used to remove low quality cells   ##### 

# hCount_RNA <- 25000
# lCount_RNA <- 1000
# hCount_ATAC <- 100000
# lCount_ATAC <- 1000
# mito_perc <- 10
# ns <- 2
# TSS.enrichment <- 1  

seurat_obj.subset <- subset(
    x = SeuratOBJ,
    subset =  (nCount_RNA > 1000 & nCount_RNA < 25000) &
        (nCount_ATAC > 1000 & nCount_ATAC < 100000)  &
        (nucleosome_signal < 2) & TSS.enrichment > 1
)

new_seurat <- paste0(base_name, '_QCed')
seurat_obj.subset$orig.ident <- new_seurat

# Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_QCed.rds'))
saveRDS(seurat_obj.subset, file = rds_name)
message('Seurat with ATAC object saved!')   

# Merge Seurat objects to get combined stats
SeuratOBJ_1.combined.GEX <- merge(SeuratOBJ, y = seurat_obj.subset,
                                  add.cell.ids = c(s_sample, new_seurat),
                                  project = "PBMC")

# Save combined Seurat object
table(SeuratOBJ_1.combined.GEX$orig.ident)
print(SeuratOBJ_1.combined.GEX)

p1_ATAC <- VlnPlot(SeuratOBJ_1.combined.GEX, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), group.by = "orig.ident") 

df_genes_per_cell <- as.data.frame(SeuratOBJ_1.combined.GEX[[]])
df_genes_per_cell %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
    #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 200, linetype=2) +
    geom_hline(yintercept = 200, linetype=2) +
    facet_wrap(~orig.ident) +
    ylab("log10(nFeature_RNA)") +
    xlab("log10(UMIs)") +
    ggtitle('UMIs/Genes by MT levels')

# Plot Genes and UMIs by density per cell 
df_genes_per_cell %>%
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +
    geom_vline(xintercept = 300) +
    ylab("Log10(UMIs)") +
    xlab("Gene-counts") +
    ggtitle("Genes density by cell") 

df_genes_per_cell %>%
    ggplot(aes(x=orig.ident, y=(nFeature_RNA), fill=orig.ident)) +
    geom_boxplot(alpha = 0.7) +
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("Log10(nFeature_RNA)") +
    xlab("") +
    ggtitle("Genes distribution by cell")



## Exercise: Merge the previous sample (pbmc3k_QC_QCed.rds) with the current FFB combined Seurat Object

## Read pre-existing Seurat for FFB 
## pre_existing_seurat <- '/Users/ccardin2/Documents/RStatClub2/processed-data/GEX_ATAC_preprocessing/'

Seurat2 <- LoadSeuratRds(here('processed-data/GEX_ATAC_preprocessing', 'pbmc3k_QC_QCed.rds'))

# Merge Seurat objects to get combined stats
table(SeuratOBJ_1.combined.GEX$orig.ident)
table(Seurat2$orig.ident)
SeuratOBJ_2.combined.GEX <- merge(SeuratOBJ_1.combined.GEX, y = Seurat2,
                                  add.cell.ids = c(s_sample, 'PBMC'),
                                  project = "PBMC-FFB")

table(SeuratOBJ_2.combined.GEX$orig.ident)

p1_ATAC <- VlnPlot(SeuratOBJ_2.combined.GEX, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), group.by = "orig.ident") 

df_genes_per_cell <- as.data.frame(SeuratOBJ_2.combined.GEX[[]])
df_genes_per_cell %>%
    ggplot(aes(x=orig.ident, y=(nFeature_RNA), fill=orig.ident)) +
    geom_boxplot(alpha = 0.7) +
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("Log10(nFeature_RNA)") +
    xlab("") +
    ggtitle("Genes distribution by cell")



############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


