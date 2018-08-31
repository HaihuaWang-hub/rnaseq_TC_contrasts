# Script to make DE Gene list for time-specific genes as id'd by DESeq2

# LOAD FILES #
deseq_tc_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_TC_genes.txt'
deseq_tc <- read.table(deseq_tc_file, header = T, stringsAsFactors = F, 
  sep = '\t')
# Arabidopsis homologs
cs_to_at_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_gene_table_v2.0.txt'
cs_to_at <- read.table(cs_to_at_file, header = T, stringsAsFactors = F, 
  sep = '\t')
# GO INFO
cs_go_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_GO_table.txt'
# Count data
count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# metadata
metadata_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt'
meta <- read.table(metadata_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #


# SET CONSTANTS #
min_floor_count_cut <- 5

# SET OUTPUT INFO #
deseq_DE_gene_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList2_DE_genes_v2.0.txt'

# LOAD PACKAGES # 
library('topGO', lib.loc = '/home/grabowsky/tools/r_packages')
library('edgeR', lib.loc = '/home/grabowsky/tools/r_packages')

#############3
deseq_tc$rank <- order(deseq_tc$adj_pval)
deseq_tc_ord <- deseq_tc[order(deseq_tc$gene), ]
deseq_DE_df <- data.frame(gene = deseq_tc_ord$gene, stringsAsFactors = F)

# Assign Arabidopsis names
gen_inds_with_at <- which(deseq_DE_df$gene %in% cs_to_at$cs_short)
at_to_incl_inds <- which(cs_to_at$cs_short %in%deseq_DE_df$gene)

deseq_DE_df$At_gene <- NA
deseq_DE_df$At_symbol <- NA
deseq_DE_df$At_description <- NA

deseq_DE_df$At_gene[gen_inds_with_at] <- cs_to_at$at_short[at_to_incl_inds]
deseq_DE_df$At_symbol[gen_inds_with_at] <- cs_to_at$at_symbol[at_to_incl_inds]
deseq_DE_df$At_description[gen_inds_with_at] <- cs_to_at$at_descr[
  at_to_incl_inds]

# Assign gene name, p-value and rank
deseq_DE_df$pval <- deseq_tc_ord$adj_pval
deseq_DE_df$rank <- deseq_tc_ord$rank

# Assign GO Terms
geneID2GO <- readMappings(file = cs_go_file)

deseq_DE_df$go_terms <- NA
for(i in seq(nrow(deseq_DE_df))){
  tmp_go_terms <- unlist(geneID2GO[deseq_DE_df$gene[i]])
  if(length(tmp_go_terms) > 0){deseq_DE_df$go_terms[i] <- paste(tmp_go_terms, 
      collapse = ';')}
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

count_inds <- which(rownames(counts_2) %in% deseq_DE_df$gene)
deseq_DE_df[, colnames(counts_2)] <- counts_2[count_inds, ]

write.table(deseq_DE_df, file = deseq_DE_gene_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

quit(save = 'no')

