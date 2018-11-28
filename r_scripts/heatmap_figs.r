# Script for generating heatmaps for the manuscript

# count_file <- args[1]
count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# note: first column of matrix should be the names of the genes 

# samp_meta_file <- args[2]
samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F,
               sep = '\t')

gene_set_1_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList1_DE_genes_v2.0.txt'
de_set1 <- read.table(gene_set_1_file, header = T, stringsAsFactors = F, 
  sep = '\t')

gene_set_2_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList2_DE_genes_v2.0.txt'
de_set2 <- read.table(gene_set_2_file, header = T, stringsAsFactors = F,
  sep = '\t')

gene_set_3_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList3_DE_genes_v2.0.txt'
de_set3 <- read.table(gene_set_3_file, header = T, stringsAsFactors = F,
  sep = '\t')

# SET VARIABLES #

# control_treatment <- args[4]
control_treatment <- 'MT5'

# SET CONSTANTS #
min_floor_count_cut <- 5
# any genes for which no libraries have a greater read count than this are
#   removed

p_cut <- 0.05
# the significance cutoff for the multiple testing-corrected p-values

# LOAD LIBRARYS #
library(DESeq2)
#library(ggplot2)
library(pheatmap)

# SET OUTPUTS #
fig_dir <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/'

geneSet1_heatmap_fig <- paste(fig_dir, 'GeneList1_heatmap.pdf', sep = '')
geneSet2_heatmap_fig <- paste(fig_dir, 'GeneList2_heatmap.pdf', sep = '')
geneSet3_heatmap_fig <- paste(fig_dir, 'GeneList3_heatmap.pdf', sep = '')

geneSet2_log2_heatmap_fig <- paste(fig_dir,
   'GeneList2_Condition_log2_heatmap.pdf', sep = '')


# SCRIPT
# adjust counts matrix
counts_1 <- counts[, -1]
rownames(counts_1) <- counts[,1]

# reorder columns so that they are in order in the heatmap
counts_1 <- counts_1[, c(7,8,9,1,2,3,4,5,6,15,16,10,11,12,13,14)]

# Re-order metadata so samples are in same order as count data
ph_order <- c()
for(i in seq(ncol(counts_1))){
  temp_ind <- which(samp_meta$SampleName == colnames(counts_1)[i])
  ph_order <- c(ph_order, temp_ind)
}
samp_meta_2 <- samp_meta[ph_order, ]

# Don't pre-normalize with DESeq2 because calculations require non-normalized
#   data

# Remove genes with low counts across samples
max_count <- apply(counts_1, 1, max)
low_rows <- which(max_count <= min_floor_count_cut)
counts_2 <- counts_1[-low_rows, ]

# format metadata
samp_meta_2$SampleName <- as.factor(samp_meta_2$SampleName)
alt_treats <- setdiff(unique(samp_meta_2$Treatment), control_treatment)
treat_order <- c(control_treatment, alt_treats)
samp_meta_2$Treatment <- factor(samp_meta_2$Treatment, levels = treat_order)
samp_meta_2$Replicate <- as.factor(samp_meta_2$Replicate)
samp_meta_2$Time <- as.factor(samp_meta_2$Time)

dds <- DESeqDataSetFromMatrix(countData = counts_2, colData = samp_meta_2,
         design = ~ Time + Treatment)
dds <- DESeq(dds)

# Variance Stabilizing Transformation
# vsd <- vst(dds, blind = FALSE)
vsd <- vst(dds)

dds_df <- as.data.frame(colData(dds)[,c('Treatment', 'Time')])
colnames(dds_df) <- c('Line', 'Time')
levels(dds_df$Line)[levels(dds_df$Line) == 'MT5'] <- 'Suneson' 
levels(dds_df$Line)[levels(dds_df$Line) == '167'] <- 'miR167OE'
levels(dds_df$Time)[levels(dds_df$Time) == '8'] <- '8 DAF'
levels(dds_df$Time)[levels(dds_df$Time) == '10'] <- '10 DAF'
levels(dds_df$Time)[levels(dds_df$Time) == '12'] <- '12 DAF'

heat_col_labs <- c(paste('miR167OE', rep(c('8 DAF', '10 DAF', '12 DAF'), 
  each = 3)), paste('Suneson', c(rep('8 DAF', times = 2), 
  rep('10 DAF', times = 2), rep('12 DAF', times = 3))))

# Gene Set 1
vsd_set1_inds <- which(names(vsd) %in% de_set1$gene)
vsd_set1 <- vsd[vsd_set1_inds]

dds_set1 <- dds[vsd_set1_inds]

vsd_set1_ord <- order(rowMeans(counts(dds_set1, normalized = T)), 
  decreasing = T)

pdf(file = geneSet1_heatmap_fig)
pheatmap(assay(vsd_set1)[vsd_set1_ord, ], cluster_rows = T, show_rownames = F,
  cluster_cols = F, annotation_col = dds_df, labels_col = heat_col_labs)
dev.off()

# Gene Set 2
vsd_set2_inds <- which(names(vsd) %in% de_set2$gene)
vsd_set2 <- vsd[vsd_set2_inds]

dds_set2 <- dds[vsd_set2_inds]

vsd_set2_ord <- order(rowMeans(counts(dds_set2, normalized = T)),
  decreasing = T)

pdf(file = geneSet2_heatmap_fig)
pheatmap(assay(vsd_set2)[vsd_set2_ord, ], cluster_rows = T, show_rownames = F,
  cluster_cols = F, annotation_col = dds_df, labels_col = heat_col_labs)
dev.off()

# Gene Set 3
vsd_set3_inds <- which(names(vsd) %in% de_set3$gene)
vsd_set3 <- vsd[vsd_set3_inds]

dds_set3 <- dds[vsd_set3_inds]

vsd_set3_ord <- order(rowMeans(counts(dds_set3, normalized = T)),
  decreasing = T)

pdf(file = geneSet3_heatmap_fig)
pheatmap(assay(vsd_set3)[vsd_set3_ord, ], cluster_rows = T, show_rownames = F,
  cluster_cols = F, annotation_col = dds_df, labels_col = heat_col_labs)
dev.off()

quit(save = 'no')


######
# Try showing log2 changes
# Is pretty ugly - can work on this if need be

ddsTC <- DESeqDataSetFromMatrix(countData = counts_2, colData = samp_meta_2,
           design = ~ Treatment + Time + Treatment:Time)
ddsTC <- DESeq(ddsTC, test = 'LRT', reduced = ~ Treatment + Time)
resTC <- results(ddsTC)

betas <- coef(ddsTC)

mat <- betas
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

pdf(file = geneSet2_log2_heatmap_fig)
pheatmap(mat[vsd_set2_inds , ], breaks = seq(from=-thr, to=thr, length = 101), 
  cluster_col = F, cluster_rows = T, show_rownames = F)
dev.off()


