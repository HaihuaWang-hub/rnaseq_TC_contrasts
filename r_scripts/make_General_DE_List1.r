# Script to make General DE gene lists from Camelina seed time course experiment
#   Get overlap of DESeq2-General, maSigPro, and splineTimeR-general

# LOAD FILES #
deseq_gen_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_general_genes.txt'
# General differences between the lines
deseq_gen <- read.table(deseq_gen_file, header = T, stringsAsFactors = F, 
  sep = '\t')

masig_tc_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_maSigPro_167_DE_genes.txt'
# Differences when accounting for time points
masig_tc <- read.table(masig_tc_file, header = T, stringsAsFactors = F, 
  sep = '\t')

splinetime_gen_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_withIntercept_genes.txt'
# All differences - overall expression level and different patterns across time
splinetime_gen <- read.table(splinetime_gen_file, header = T, 
  stringsAsFactors = F, sep = '\t')
## File with A.thaliana homologs
cs_to_at_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_gene_table_v2.0.txt'
cs_to_at <- read.table(cs_to_at_file, header = T, stringsAsFactors = F, 
  sep = '\t')
# GO info
cs_go_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_GO_table.txt'
# Count data
count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# Metadata
metadata_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt'
meta <- read.table(metadata_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #


# SET CONSTANTS #
min_floor_count_cut <- 5

# SET OUTPUT INFO #
gen_DE_gene_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList1_DE_genes_v2.0.txt'

# LOAD PACKAGES #
library('topGO', lib.loc = '/home/grabowsky/tools/r_packages')
library('edgeR', lib.loc = '/home/grabowsky/tools/r_packages')

####################3
# add rank columns to result files
deseq_gen$rank <- order(deseq_gen$padj)
masig_tc$rank <- order(masig_tc$p_val)
splinetime_gen$rank <- order(splinetime_gen$adj_pval)

# Generate "General" Gene list
gen_names <- intersect(deseq_gen$gene, intersect(masig_tc$genes, splinetime_gen$gene))

# Make file with "general" genes
gen_names_2 <- sort(gen_names)
deseq_gen_ord <- deseq_gen[order(deseq_gen$gene),]
masig_tc_ord <- masig_tc[order(masig_tc$genes),]
splinetime_gen_ord <- splinetime_gen[order(splinetime_gen$gene), ]

# p-values are all multiple-testing corrected
gen_DE_df <- data.frame(gene = gen_names_2, stringsAsFactors = F)

# ASSIGN ARABIDOPSIS NAMES
gen_inds_with_at <- which(gen_DE_df$gene %in% cs_to_at$cs_short)
at_to_incl_inds <- which(cs_to_at$cs_short %in%gen_DE_df$gene)

gen_DE_df$At_gene <- NA
gen_DE_df$At_symbol <- NA
gen_DE_df$At_description <- NA

gen_DE_df$At_gene[gen_inds_with_at] <- cs_to_at$at_short[at_to_incl_inds]
gen_DE_df$At_symbol[gen_inds_with_at] <- cs_to_at$at_symbol[at_to_incl_inds]
gen_DE_df$At_description[gen_inds_with_at] <- cs_to_at$at_descr[at_to_incl_inds]

# ASSIGN GENE NAME, P-VALUES AND RANKS FOR EACH APPROACH
gen_DE_df$deseq_pval <- deseq_gen_ord$padj[which(deseq_gen_ord$gene %in% 
  gen_DE_df$gene)]
gen_DE_df$deseq_rank <- deseq_gen_ord$rank[which(deseq_gen_ord$gene %in% 
  gen_DE_df$gene)]
gen_DE_df$masig_pval <- masig_tc_ord$p_val[which(masig_tc_ord$genes %in% 
  gen_DE_df$gene)]
gen_DE_df$masig_rank <- masig_tc_ord$rank[which(masig_tc_ord$genes %in% 
  gen_DE_df$gene)]
gen_DE_df$spline_pval <- splinetime_gen_ord$adj_pval[which(
  splinetime_gen_ord$gene %in% gen_DE_df$gene)]
gen_DE_df$spline_rank <- splinetime_gen_ord$rank[which(
  splinetime_gen_ord$gene %in% gen_DE_df$gene)]

# NOTES:
# lowest missing rank for deseq is 73; missing 73, 84, 94, and 99 from top 100
# for maSigPro, only has 4 and 7 from top 10, missing 53 of top 100
# lowest missing rank for splineTimeR is 9, missing 9 of top 100

# ASSIGN GO TERMS
geneID2GO <- readMappings(file = cs_go_file)

gen_DE_df$go_terms <- NA
for(i in seq(nrow(gen_DE_df))){
  tmp_go_terms <- unlist(geneID2GO[gen_DE_df$gene[i]])
  if(length(tmp_go_terms) > 0){gen_DE_df$go_terms[i] <- paste(tmp_go_terms, collapse = ';')}
}

# ASSIGN EXPRESSION LEVELS
counts_1 <- counts[, -1]
rownames(counts_1) <- counts[,1]

norm_factors <- calcNormFactors(counts_1, method = 'TMM')
counts_2 <- mapply(`/`, counts_1, norm_factors)
rownames(counts_2) <- rownames(counts_1)

for(i in seq(ncol(counts_2))){
  meta_ind <- which(meta$SampleName == colnames(counts_2)[i])
  colnames(counts_2)[i] <- paste(meta$sampleName[meta_ind], meta$SampleName[meta_ind], sep = '_')
}

count_inds <- which(rownames(counts_2) %in% gen_DE_df$gene)
gen_DE_df[, colnames(counts_2)] <- counts_2[count_inds, ]

write.table(gen_DE_df, file = gen_DE_gene_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)


quit(save = 'no')

