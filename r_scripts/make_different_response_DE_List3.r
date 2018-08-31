# Script to make List 3: Different response across time

# LOAD FILES #
splinetime_tc_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_noIntercept_genes.txt'
splinetime_tc <- read.table(splinetime_tc_file, header = T, stringsAsFactors = F, sep = '\t')
# Arabidopsis homologs
cs_to_at_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_gene_table_v2.0.txt'
cs_to_at <- read.table(cs_to_at_file, header = T, stringsAsFactors = F, 
  sep = '\t')
# GO Info
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
splinetime_DE_gene_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList3_DE_genes_v2.0.txt'

# LOAD LIBRARIES #
library('topGO', lib.loc = '/home/grabowsky/tools/r_packages')
library('edgeR', lib.loc = '/home/grabowsky/tools/r_packages')

#####################3
splinetime_tc$rank <- order(splinetime_tc$adj_pval)
splinetime_tc_ord <- splinetime_tc[order(splinetime_tc$gene), ]
splinetime_DE_df <- data.frame(gene = splinetime_tc_ord$gene, 
  stringsAsFactors = F)

# Assign Arabidopsis Names
gen_inds_with_at <- which(splinetime_DE_df$gene %in% cs_to_at$cs_short)
at_to_incl_inds <- which(cs_to_at$cs_short %in%splinetime_DE_df$gene)

splinetime_DE_df$At_gene <- NA
splinetime_DE_df$At_symbol <- NA
splinetime_DE_df$At_description <- NA

splinetime_DE_df$At_gene[gen_inds_with_at] <- cs_to_at$at_short[at_to_incl_inds]
splinetime_DE_df$At_symbol[gen_inds_with_at] <- cs_to_at$at_symbol[at_to_incl_inds]
splinetime_DE_df$At_description[gen_inds_with_at] <- cs_to_at$at_descr[at_to_incl_inds]
# Assign p-values, and ranks
splinetime_DE_df$pval <- splinetime_tc_ord$adj_pval
splinetime_DE_df$rank <- splinetime_tc_ord$rank

# Assign GO Terms
geneID2GO <- readMappings(file = cs_go_file)

splinetime_DE_df$go_terms <- NA
for(i in seq(nrow(splinetime_DE_df))){
  tmp_go_terms <- unlist(geneID2GO[splinetime_DE_df$gene[i]])
  if(length(tmp_go_terms) > 0){splinetime_DE_df$go_terms[i] <- paste(tmp_go_terms, collapse = ';')}
}

# Assign Expression Levels
counts_1 <- counts[, -1]
rownames(counts_1) <- counts[,1]

norm_factors <- calcNormFactors(counts_1, method = 'TMM')
counts_2 <- mapply(`/`, counts_1, norm_factors)
rownames(counts_2) <- rownames(counts_1)

for(i in seq(ncol(counts_2))){
  meta_ind <- which(meta$SampleName == colnames(counts_2)[i])
  colnames(counts_2)[i] <- paste(meta$sampleName[meta_ind], 
    meta$SampleName[meta_ind], sep = '_')
}

count_inds <- which(rownames(counts_2) %in% splinetime_DE_df$gene)
splinetime_DE_df[, colnames(counts_2)] <- counts_2[count_inds, ]

# Write table
write.table(splinetime_DE_df, file = splinetime_DE_gene_out_file, quote = F, 
  sep = '\t', row.names = F, col.names = T)

quit(save = 'no')

