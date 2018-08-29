# Script for identifying differentially expressed genes in a time course
#  experiment of two conditions (Control vs Treatment) using the splineTimeR 
#  package

# Arguments:
# [1]: (character) filename of read-count expression matrix; rows = genes, 
#        columns = read counts per library/sample; first column = gene names;
#        column names must have sample/library name
# [2]: (character) filename of metadata for sequencing libraries;
#        Must contain:
#        i) column named <SampleName> with names of technical samples/libraries
#             that corresponds to colunn names of read-count matrix
#        ii) column named <Treatment> which indicates which overall biological
#              sample, line, or treatment the library comes from 
#              ex: Control, Salt, 167 
# [3]: (character) prefix for the file names that will be returned/saved
# [4]: (character) the name of the Control line or treatment as written
#        in the <Treatment> column of the sample metadata
#################

args = commandArgs(trailingOnly = TRUE)

# INPUT FILES #
count_file <- args[1] 
###count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# NOTE: first column of matrix should be the names of the genes 

samp_meta_file <- args[2]
###samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F,
               sep = '\t')

# SET VARIABLES #
out_pre <- args[3]
###out_pre <- 'test'

control_treatment <- args[4]
###control_treatment <- 'MT5'

# SET CONSTANTS #
min_floor_count_cut <- 5
# any genes for which no libraries have a greater read count than this are
#   removed

# LOAD LIBRARYS #
library('edgeR', lib.loc = '/home/grabowsky/tools/r_packages')
library('Biobase', lib.loc = '/home/grabowsky/tools/r_packages')
# following libraries are for splineTimerR
library(FIs, lib.loc = '/home/grabowsky/tools/r_packages')
library(fdrtool, lib.loc = '/home/grabowsky/tools/r_packages')
library(longitudinal, lib.loc = '/home/grabowsky/tools/r_packages')
library(GeneNet, lib.loc = '/home/grabowsky/tools/r_packages')
library(splineTimeR, lib.loc = '/home/grabowsky/tools/r_packages')

#####################
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

# Normalize read counts using edgeR function
norm_factors <- calcNormFactors(counts_1, method = 'TMM')
counts_2 <- mapply(`/`, counts_1, norm_factors)
rownames(counts_2) <- rownames(counts_1)

# Remove genes with low counts across samples
max_count <- apply(counts_2, 1, max)
low_rows <- which(max_count <= min_floor_count_cut)
counts_3 <- counts_2[-low_rows, ]
# REMOVE FOLLOWING LINE WHEN DONE TESTING!!!!
#counts_3 <- counts_3[c(1:500), ]

# format metadata
samp_meta_2$SampleName <- as.factor(samp_meta_2$SampleName)
samp_meta_2$Treatment <- as.factor(samp_meta_2$Treatment)
samp_meta_2$Replicate <- as.factor(samp_meta_2$Replicate)
samp_meta_2$Time <- as.numeric(samp_meta_2$Time)
rownames(samp_meta_2) <- samp_meta_2$SampleName

meta_annotDF <- new('AnnotatedDataFrame', samp_meta_2)

# Generate ExpressionSet Object required for splineTimeR
expr_set <- ExpressionSet(assayData = counts_3, phenoData = meta_annotDF)

# Make list of DE results
# Element 1: DE genes including those with baseline differential expression;
# Element 2: DE genes that only show different response to time
de_list <- list()
de_list[['withIntercept']] <- splineDiffExprs(eSetObject = expr_set, df = 2, 
  cutoff.adj.pVal = 0.05, reference = control_treatment, intercept = T)
#
de_list[['noIntercept']] <- splineDiffExprs(eSetObject = expr_set, df = 2, 
  cutoff.adj.pVal = 0.05, reference = control_treatment, intercept = F)

# adjust gene names to 'character' and adjust column names for output files
for(i in seq(length(de_list))){
  de_list[[i]]$row_IDs <- as.character(de_list[[i]]$row_IDs)
  colnames(de_list[[i]])[which(colnames(de_list[[i]]) == 'row_IDs')] <-
  'gene'
  colnames(de_list[[i]])[which(colnames(de_list[[i]]) == 'P.Value')] <-
  'p_val'
  colnames(de_list[[i]])[which(colnames(de_list[[i]]) == 'adj.P.Val')] <-
  'adj_pval'
}

# Output files 
# for each gene set: 1) gene names and p-values; 2) gene names and spline
#  values; 3) full matrix
for(j in seq(length(de_list))){
  data_name <- names(de_list)[j]
  gene_file_name <- paste( paste(out_pre, 'splineTimeR', data_name, 
    'genes', sep = '_'), '.txt', sep = '')
  write.table(de_list[[j]][, c('gene', 'p_val', 'adj_pval')], 
    file = gene_file_name, quote = F, sep = '\t', row.names = F, col.names = T)
  #
  paramval_file_name <- paste( paste(out_pre, 'splineTimeR', data_name, 
   'paramVals', sep = '_'), '.txt', sep = '')
  p_rm_cols <- c('AveExpr', 'F', 'p_val', 'adj_pval')
  p_rm_inds <- c()
  for(k in p_rm_cols){
    tmp_ind <- which(colnames(de_list[[j]]) == k)
    p_rm_inds <- c(p_rm_inds, tmp_ind)
  }
  write.table(de_list[[j]][, p_rm_inds], file = paramval_file_name, quote = F, 
    sep = '\t', row.names = F, col.names = T)
  #
  full_file_name <- paste( paste(out_pre, 'splineTimeR', data_name, 
    'full_mat', sep = '_'), '.txt', sep = '')
  write.table(de_list[[j]], file = full_file_name, quote = F, sep = '\t', 
    row.names = F, col.names = T)
}

quit(save = 'no')
