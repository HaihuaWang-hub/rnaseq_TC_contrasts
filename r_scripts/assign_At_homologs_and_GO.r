# Script for assigning A.thaliana homologs and GO terms to Camelina genes using
#   the results from BLASTP allignments using <diamond> and GO terms from
#   Arabidopsis

# LOAD FILES #
cam_hit_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/cam_arab_matches.m8'
cam_hits <- read.table(cam_hit_file, header = F, stringsAsFactors = F)

at_go_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/ATH_GO_GOSLIM.txt'
at_go <- scan(at_go_file, what = 'character', sep = '$', quote = '')
at_go_1 <- strsplit(at_go, split = '\t')

# SET VARIABLES #

# SET CONSTANTS #

# SET OUTPUT NAMES #
# For full table with more info
tot_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_gene_table.txt'

# For table with just Camelina gene and GO terms
cs_go_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_GO_table.txt'

# LOAD PACKAGES #


######################
# Assign homologs
hit_col_names <- c('cs_gene', 'at_gene', 'cs_len', 'at_length', 'eval', 
  'bit', 'match_len', 'per_ident')
colnames(cam_hits) <- hit_col_names

cam_sing_hits <- cam_hits[-which(duplicated(cam_hits$cs_gene)), ]
cam_sing_hits <- cam_sing_hits[order(cam_sing_hits$cs_gene), ]
cam_sing_hits$n_hits <- table(cam_hits$cs_gene)

cam_sing_hits$go_terms <- NA
tmp_at_gene <- strsplit(cam_sing_hits$at_gene, split = '.', fixed = T)
cam_sing_hits$at_short <- unlist(lapply(tmp_at_gene, function(x) x[1]))

tmp_cs_gene <- strsplit(cam_sing_hits$cs_gene, split = '.', fixed = T)
cam_sing_hits$cs_short <- unlist(lapply(tmp_cs_gene, function(x) x[1]))

# Assign GO Terms
## Process Arabidopsis GO Terms
at_gene_vec <- unlist(lapply(at_go_1, function(x) x[[1]]))
at_go_vec <- unlist(lapply(at_go_1, function(x) 
  x[grep('GO:', x, fixed = T)[1]]))

at_gene_inds <- grep('^AT[0-9]', at_gene_vec, fixed = F)
at_gene_vec_1 <- at_gene_vec[at_gene_inds]
at_go_vec_1 <- at_go_vec[at_gene_inds] 

gene_order <- order(at_gene_vec_1)
at_gene_vec_2 <- at_gene_vec_1[gene_order]
at_go_vec_2 <- at_go_vec_1[gene_order]
at_gene_go <- tapply(at_go_vec_2, at_gene_vec_2, function(x) 
  paste(unique(x), collapse = ', '))
at_go_df <- data.frame(at_gene = names(at_gene_go), go_terms = at_gene_go, 
  stringsAsFactors = F)

at_in_vec <- which(at_go_df$at_gene %in% cam_sing_hits$at_short)
at_go_df_2 <- at_go_df[at_in_vec, ]

## Match Arabidopsis GO terms with Camelina genes
for(i in seq(nrow(at_go_df_2))){
  tmp_inds <- which(cam_sing_hits$at_short == at_go_df_2$at_gene[i])
  cam_sing_hits$go_terms[tmp_inds] <- at_go_df_2$go_terms[i]
  print(i)
}

max_reps <- max(table(cam_sing_hits$cs_short)) 
gene_tab_names <- names(table(cam_sing_hits$cs_short))
rm_inds <- c()

for(i in c(2:max_reps)){
  rep_genes <- gene_tab_names[which(table(cam_sing_hits$cs_short) == i)]
  for(g in rep_genes){
    tmp_inds <- which(cam_sing_hits$cs_short == g)
    top_ind <- which.max(cam_sing_hits$bit[tmp_inds])
    tmp_rm_inds <- tmp_inds[-top_ind]
    rm_inds <- c(rm_inds, tmp_rm_inds)
  }
}

cam_sing_hits_filt <- cam_sing_hits[-rm_inds, ]

write.table(cam_sing_hits_filt, file = tot_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

cs_go_df <- cam_sing_hits_filt[, c('cs_short', 'go_terms')]
write.table(cs_go_df, file = cs_go_out_file, quote = F, sep = '\t', row.names = F, col.names = F)

quit(save = 'no')

