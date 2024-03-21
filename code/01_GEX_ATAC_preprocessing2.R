########################################################################
##
## RstatClub: Multiome Analysis with Seurat and Signac #2
##
## Input: Truncated H5, meta-data.cvs and fragments.tvs
## Output: rds Seurat objects with ATAC assay, one with QCs and other w/o QCs
##.        GEX and ATAC plots
##
## NOTES: 
## For a ~10k cells it is recommended ~40G free mem to process the TSS() ATAC score.  
## Without storing the base-resolution matrix of integration counts at each site you can use less memory
## For slurm env: $srun --pty --mem=40GB --x11 bash
##
########################################################################

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


########################    Initials.  ########################  

## Link to dataset: https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0 
## Flash-Frozen Human Healthy Brain Tissue (3k)

sample_sample <- 'FFB_Healthy'
#s_tissue <- 'human'

message('Processing sample: ',sample_sample, ' from ', s_tissue, ' brain tissue.')

base_path <- paste0('processed-data/cellranger-arc')    # Read H5 file
meta_path <- here(base_path, "human_brain_3k_per_barcode_metrics.csv")
bc_mtx_path <- here(base_path, "human_brain_3k_filtered_feature_bc_matrix.h5" )
atac_path <- here(base_path,"human_brain_3k_atac_fragments.tsv.gz")

    
# filtered bc mtx
mtx <- Read10X_h5(bc_mtx_path)             # Returns a sparse matrix with rows and columns labeled
rna_counts <- mtx$`Gene Expression`

metadata <- read.csv(file = meta_path, header = TRUE, row.names = 1)
# subset specific fields in the meta data df
meta_tmp = c('atac_peak_region_fragments','atac_fragments')
meta = metadata[meta_tmp]
print(head(meta, n = 3))

SeuratOBJ <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = sample_sample,
    meta.data = meta)

## Data verification
SeuratOBJ
str(SeuratOBJ)
head(SeuratOBJ, n=3)
Idents(SeuratOBJ)

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(sample_sample,'.rds'))
saveRDS(SeuratOBJ, file = rds_name)

## Read pre-existing Seurat
## pre_existing_seurat <- '/Users/ccardin2/Documents/RStatClub2/processed-data/GEX_ATAC_preprocessing/'
#Seurat2 <- LoadSeuratRds(here('processed-data/GEX_ATAC_preprocessing', paste0(sample_sample,'.rds')))

## Calculate MITO levels
SeuratOBJ$log10GenesPerUMI <- log10(SeuratOBJ$nFeature_RNA) / log10(SeuratOBJ$nCount_RNA)

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
p1_GEX <- VlnPlot(SeuratOBJ, features = sfeature, group.by = 'orig.ident')
png_file <- paste0(base_name,'_GEX_QCs.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(p1_GEX, filename = png_name, height = 4, width = 4)

## Plot Genes and UMIs by density per cell 

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




########  Create ATAC assay and attach it to Seurat object     ##### 

message('Processing ATAC for sample ', levels(SeuratOBJ$`orig.ident`[1]))


## Create gene annotations for hg38 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# show(annotations) / # View(head(annotations,n=5)) /# names(genomeStyles('Homo_sapiens'))

## Format annotation 
seqlevelsStyle(annotations) <- "UCSC"
length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) 
genome(annotations) <- "hg38"
# names(mcols(annotations))

## Base name to save plots
base_name <- paste0(levels(SeuratOBJ$`orig.ident`[1]),'_QC')

## Pull peak fragment information and 
atac_counts <- mtx$Peaks

# Convert a genomic coordinate string to a GRanges object and attach to ATAC counts
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))      
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Create chromatin assay
chrom_assay <- CreateChromatinAssay(counts = atac_counts, 
                                    sep = c(":", "-"), 
                                    fragments = atac_path, 
                                    annotation = annotations)


SeuratOBJ[["ATAC"]] <- chrom_assay

## Annotations of the object are set
Annotation(SeuratOBJ[["ATAC"]]) <- annotations

head(SeuratOBJ, n = 3)

message('Chromatin assay attached to seurat object successfully')



## Measure ATAC quality matrics:  NFR and TSS

DefaultAssay(SeuratOBJ) <- "ATAC"

## Calculate the strength of the nucleosome signal per cell

SeuratOBJ <- NucleosomeSignal(SeuratOBJ)
SeuratOBJ$nucleosome_group <- ifelse(SeuratOBJ$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p1_NS <- FragmentHistogram(object = SeuratOBJ, group.by = 'nucleosome_group')

png_file <- paste0(base_name,'_NS.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(p1_NS, filename = png_name, height = 4, width = 4)

## Calculate the "Transcription Start Site (TSS)" enrichment score

tryCatch( {
    # Important NOTES: https://github.com/stuart-lab/signac/issues/374 
    # 1. Some times jumps an issue of ** memory allocation ** when "fast=FALSE", FALSE argument is mandatory to compute the TSS enrichment scores and visualize with TSSPlot(). You need more memory
    
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
    
    ## Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$blacklist_fraction <- FractionCountsInRegion(
        object = SeuratOBJ,
        assay = 'ATAC',
        regions = blacklist_hg38)
    
    ## Add blacklist ratio and fraction of reads in peaks
    SeuratOBJ$pct_reads_in_peaks <- SeuratOBJ$atac_peak_region_fragments / SeuratOBJ$atac_fragments * 100
    SeuratOBJ$blacklist_ratio <- SeuratOBJ$blacklist_fraction / SeuratOBJ$atac_peak_region_fragments
    #        Error in `x[[i, drop = TRUE]]`:
    #       ! 'blacklist_fraction' not found in this Seurat object
    
    # Plot Peaks in black ratio 
    p1_BlackR <-  VlnPlot(SeuratOBJ, features = c("pct_reads_in_peaks","blacklist_ratio"), group.by = "orig.ident", ncol = 2)  
    png_file_BlackR <- paste0(base_name,'_reads_in_peaks.png')
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_BlackR)
    ggsave(p1_BlackR, filename = png_name, height = 4, width = 5)
    
    ## Plot ATAC main feature scores
    p1_ATAC <- VlnPlot(SeuratOBJ, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), group.by = "orig.ident") 
    png_file_ATAC <- paste0(base_name,'_ATAC_QCs.png')
    png_name <- here('plots/GEX_ATAC_preprocessing', png_file_ATAC)
    ggsave(p1_ATAC, filename = png_name, height = 4, width = 7)
    
    message('ATAC QCs plots saved!')  
    

# Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_GEX_ATAC.rds'))
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC object saved!')   



########  Standard thresholds to remove low quality cells   ##### 

## Apply arbitrary scores as example

SeuratOBJ_QCed <- subset(x = SeuratOBJ,
    subset =  (nCount_RNA > 1000 & nCount_RNA < 10000) &
        (nCount_ATAC > 1000 & nCount_ATAC < 100000)  &
        (nucleosome_signal < 2) & TSS.enrichment > 1)

new_seurat <- paste0(base_name, '_QCed')
SeuratOBJ_QCed$orig.ident <- new_seurat

## Save Seurat QCed data

rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_QCed.rds'))
saveRDS(SeuratOBJ_QCed, file = rds_name)


message('Seurat QCed saved!')   






# # Merge Seurat objects to get combined stats
# SeuratOBJ_1.combined.GEX <- merge(SeuratOBJ, y = seurat_obj.subset,
#                                   add.cell.ids = c(sample_sample, new_seurat),
#                                   project = "PBMC")
# 
# # Save combined Seurat object
# table(SeuratOBJ_1.combined.GEX$orig.ident)
# print(SeuratOBJ_1.combined.GEX)
# 
# p1_ATAC <- VlnPlot(SeuratOBJ_1.combined.GEX, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), group.by = "orig.ident") 
# 
# df_genes_per_cell <- as.data.frame(SeuratOBJ_1.combined.GEX[[]])
# df_genes_per_cell %>%
#     ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
#     #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
#     geom_point() +
#     scale_colour_gradient(low = "gray90", high = "black") +
#     stat_smooth(method=lm) +
#     scale_x_log10() +
#     scale_y_log10() +
#     theme_classic() +
#     geom_vline(xintercept = 200, linetype=2) +
#     geom_hline(yintercept = 200, linetype=2) +
#     facet_wrap(~orig.ident) +
#     ylab("log10(nFeature_RNA)") +
#     xlab("log10(UMIs)") +
#     ggtitle('UMIs/Genes by MT levels')
# 
# # Plot Genes and UMIs by density per cell 
# df_genes_per_cell %>%
#     ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
#     geom_density(alpha = 0.2) +
#     scale_x_log10() +
#     theme_classic() +
#     theme(plot.title = element_text(hjust=0.5)) +
#     geom_vline(xintercept = 300) +
#     ylab("Log10(UMIs)") +
#     xlab("Gene-counts") +
#     ggtitle("Genes density by cell") 
# 
# df_genes_per_cell %>%
#     ggplot(aes(x=orig.ident, y=(nFeature_RNA), fill=orig.ident)) +
#     geom_boxplot(alpha = 0.7) +
#     theme_classic() +
#     theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
#     theme(plot.title = element_text(hjust=0.5)) +
#     ylab("Log10(nFeature_RNA)") +
#     xlab("") +
#     ggtitle("Genes distribution by cell")




############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


