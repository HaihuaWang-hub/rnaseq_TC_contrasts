# Script to add the Arabidopsis common name info for the homologs assigned
#  to the Camelina genes

# LOAD FILES #
tot_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_gene_table.txt'
tot_tab <- read.table(tot_file, header = T, stringsAsFactors = F, sep = '\t')

alias_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/gene_aliases_20130831.txt'
gene_alias <- scan(alias_file, what = 'character', sep = '\n')

# SET VARIABLES #

# SET CONSTANTS #

# SET OUTPUT NAMES #
# Full full table
tot_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_gene_table_v2.0.txt'

# Table with less info
pert_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_to_AT_info_simple_table.txt'

# LOAD PACKAGES #


###############
tot_ord <- tot_tab[order(tot_tab$at_short), ]

alias_list <- strsplit(gene_alias, split = '\t')
alias_length <- lapply(alias_list, length)
table(unlist(alias_length))

alias_at_genes <- unlist(lapply(alias_list, function(x) x[1]))[-1]
alias_symbol <- unlist(lapply(alias_list, function(x) x[2]))[-1]
alias_fullname <- rep(NA, times = length(alias_list))
alias_fullname[which(alias_length == 3)] <- unlist(
  lapply(alias_list[which(alias_length == 3)], function(x) x[3]))
alias_fullname <- alias_fullname[-1]

alias_df <- data.frame(
  gene = alias_at_genes[which(duplicated(alias_at_genes) == F)], 
  symbol = alias_symbol[which(duplicated(alias_at_genes) == F)], 
  fullname = alias_fullname[which(duplicated(alias_at_genes) == F)], 
  stringsAsFactors = F)

dup_genes_df <- data.frame(
  gene = alias_at_genes[which(duplicated(alias_at_genes))], 
  symbol = alias_symbol[which(duplicated(alias_at_genes))], 
  fullname = alias_fullname[which(duplicated(alias_at_genes))], 
  stringsAsFactors = F)

for(i in seq(nrow(dup_genes_df))){
  df_ind <- which(alias_df$gene == dup_genes_df$gene[i])
  alias_df$symbol[df_ind] <- paste(alias_df$symbol[df_ind], 
                               dup_genes_df$symbol[i], sep = ';')
  if(is.na(dup_genes_df$fullname[i]) == F){
    tmp_fullname <- alias_df$fullname[df_ind]
    new_fullname <- paste(tmp_fullname, dup_genes_df$fullname[i], sep = ';')
    if(is.na(tmp_fullname)){new_fullname <- dup_genes_df$fullname[i]}
    alias_df$fullname[df_ind] <- new_fullname
  }
  print(i)
}

# from here, match with the Camelina dataframe

tot_in_alias_inds <- which(tot_ord$at_short %in% alias_df$gene)
alias_in_tot_inds <- which(alias_df$gene %in% tot_ord$at_short)

tot_ord$at_symbol <- NA
tot_ord$at_descr <- NA

for(i in alias_in_tot_inds){
  tmp_tot_inds <- which(tot_ord$at_short == alias_df$gene[i])
  tot_ord$at_symbol[tmp_tot_inds] <- alias_df$symbol[i]
  tot_ord$at_descr[tmp_tot_inds] <- alias_df$fullname[i]
  print(i)
}

tot_reord <- tot_ord[order(tot_ord$cs_short), ]
tot_reord$at_descr <- gsub("\"", '', tot_reord$at_descr)
tot_reord$at_descr <- gsub('5\'', '5-prime', tot_reord$at_descr)
tot_reord$at_descr <- gsub('3\'', '3-prime', tot_reord$at_descr)
tot_reord$at_descr <- gsub('C\'', 'C-prime', tot_reord$at_descr)
tot_reord$at_descr <- gsub('2\'', '2-prime', tot_reord$at_descr)
tot_reord$at_descr <- gsub('\'', '', tot_reord$at_descr)
tot_reord$at_descr <- gsub('#', '', tot_reord$at_descr)

tot_reord$at_symbol <- gsub('\'', '-prime', tot_reord$at_symbol)
tot_reord$at_symbol <- gsub('#', '', tot_reord$at_symbol)

write.table(tot_reord, file = tot_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

# Generate file that has Camelina gene (in format used in analyzed data), At gene, At symbol, and At description
pert_info <- tot_reord[, c('cs_short', 'at_short', 'at_symbol', 'at_descr')]
write.table(pert_info, file = pert_out_file, quote = F, sep = '\t', row.names = F, col.names = T)

quite(save = 'no')

