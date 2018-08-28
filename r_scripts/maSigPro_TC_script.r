# Script for identifying differentially expressed genes between conditions
#  across a time course experiment using the maSigPro R package

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
library('maSigPro', lib.loc = '/home/grabowsky/tools/r_packages')

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

# generate edisign object
samp_meta_2$Condition <- as.numeric( factor( paste( samp_meta_2$Treatment, 
                           samp_meta_2$Time, sep = '_')))
edes <- data.frame(Time = samp_meta_2$Time, Replicate = samp_meta_2$Condition, 
          stringsAsFactors = F)
rownames(edes) <- samp_meta_2$SampleName

edes$Control <- as.numeric(samp_meta_2$Treatment == control_treatment)
alt_treat_vec <- setdiff(unique(samp_meta_2$Treatment), control_treatment)

for(t in alt_treat_vec){
  treat_name <- paste('treament', t, sep = '_')
  edes[,treat_name] <- as.numeric(samp_meta_2$Treatment == t)
}

design <- make.design.matrix(edes, degree = 2)

# Normalize read counts using edgeR function
norm_factors <- calcNormFactors(counts_1, method = 'TMM')
counts_2 <- mapply(`/`, counts_1, norm_factors)
rownames(counts_2) <- rownames(counts_1)

# Remove genes with low counts across samples
max_count <- apply(counts_2, 1, max)
low_rows <- which(max_count <= min_floor_count_cut)
counts_3 <- counts_2[-low_rows, ]

# estimate the GLM theta using edgeR function
est_theta <- 1/(estimateGLMCommonDisp(counts_3))

# run functions to find significant genes
fit <- p.vector(counts_3, design, Q = 0.05, MT.adjust = 'BH', counts = T, 
         min.obs = 1, theta = est_theta)
tstep <- T.fit(fit, step.method = 'backward', alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.7, vars = 'groups')

# Generate output files for the differnt treatments; no files generated for
#   "Control" treatment

file_types_names <- c('DE_genes', 'DE_betas', 'DE_full_mat')
for(i in seq(length(alt_treat_vec))){
  #i <- 1
  out_list <- list()
  beta_tab <- sigs$sig.genes[[(i+1)]]$coefficients
  pval_tab <- sigs$sig.genes[[(i+1)]]$sig.pvalues
  # data.frame with gene names and overall p-value
  out_list[['DE_genes']] <- data.frame(genes = rownames(pval_tab), 
                              p_val = pval_tab$'p-value', stringsAsFactors = F)
  # data.frame with betas for all variables; can try to use for clustering
  #   by behaviour across conditions
  out_list[['DE_betas']] <- data.frame(genes = rownames(beta_tab), beta_tab, 
                              stringsAsFactors = F)
  # data.frame with overall p-value, R^2, betas for each variable, and
  #   p-values for each beta
  out_list[['DE_full_mat']] <- data.frame(genes = rownames(pval_tab), 
                                 pval_tab, beta_tab, stringsAsFactors = F)
  for(df in seq(length(out_list))){
    out_name <- paste( paste(out_pre, alt_treat_vec[i], 
                  names(out_list)[df], sep = '_'), '.txt', sep = '')
    write.table(out_list[[df]], file = out_name, quote = F, sep = '\t', 
      row.names = F, col.names = T)
  }
}

quit(save = 'no')
