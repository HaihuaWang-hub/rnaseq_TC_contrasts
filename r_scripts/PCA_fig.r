# Script for generating PCA plot for Camelina manuscript
# Take-homes: 
#  Tried using both VST and rlog transformations - they both produce
#    almost the exact same plots.
#  I'll provide the plot based on VST since that's faster to re-generate, if 
#    need be.

# count_file <- args[1]
count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# note: first column of matrix should be the names of the genes 

# samp_meta_file <- args[2]
samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F,
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
library(ggplot2)

# SET OUTPUTS #
fig_dir <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/'

vst_pca_file <- paste(fig_dir, 'PCA_vst.pdf', sep = '')
rlog_pca_file <- paste(fig_dir, 'PCA_rlog.pdf', sep = '')

# SCRIPT
# adjust counts matrix
counts_1 <- counts[, -1]
rownames(counts_1) <- counts[,1]

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

vsd_pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Time"), returnData = T)
v_percVar <- round(100 * attr(vsd_pcaData, "percentVar"))
v_ggPCA <- ggplot(vsd_pcaData, aes(PC1, PC2, color = Treatment, shape = Time)) +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', v_percVar[1], '% variance')) +
  ylab(paste0('PC2: ', v_percVar[2], '% variance'))
#  coord_fixed()

# Adjust labels and names for better description in plot
v_ggPCA$labels$colour <- 'Line'

levels(v_ggPCA$data$Treatment)[
  levels(v_ggPCA$data$Treatment)=='MT5'] <- 'Suneson'
levels(v_ggPCA$data$Treatment)[
  levels(v_ggPCA$data$Treatment)=='167'] <- 'miR167OE'

levels(v_ggPCA$data$Time)[
  levels(v_ggPCA$data$Time)=='8'] <- '8 DAF'
levels(v_ggPCA$data$Time)[
  levels(v_ggPCA$data$Time)=='10'] <- '10 DAF'
levels(v_ggPCA$data$Time)[
  levels(v_ggPCA$data$Time)=='12'] <- '12 DAF'


pdf(file = vst_pca_file, height = 4, width = 5.5)
v_ggPCA
dev.off()

###########
# regularized log transformation
# rld <- rlog(dds, blind = F)
rld <- rlog(dds)

r_pcaData <- plotPCA(rld, intgroup=c('Treatment', 'Time'), returnData = T)
r_percVar <- round(100 * attr(r_pcaData, "percentVar"))
r_ggPCA <- ggplot(r_pcaData, aes(PC1, PC2, color = Treatment, shape = Time)) +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', r_percVar[1], '% variance')) +
  ylab(paste0('PC2: ', r_percVar[2], '% variance'))
#  coord_fixed()

# Adjust labels and names for better description in plot
r_ggPCA$labels$colour <- 'Line'

levels(r_ggPCA$data$Treatment)[
  levels(r_ggPCA$data$Treatment)=='MT5'] <- 'Suneson'
levels(r_ggPCA$data$Treatment)[
  levels(r_ggPCA$data$Treatment)=='167'] <- 'miR167OE'

levels(r_ggPCA$data$Time)[
  levels(r_ggPCA$data$Time)=='8'] <- '8 DAF'
levels(r_ggPCA$data$Time)[
  levels(r_ggPCA$data$Time)=='10'] <- '10 DAF'
levels(r_ggPCA$data$Time)[
  levels(r_ggPCA$data$Time)=='12'] <- '12 DAF'

pdf(file = rlog_pca_file, height = 4, width = 5.5)
r_ggPCA
dev.off()


quit(save = 'no')

