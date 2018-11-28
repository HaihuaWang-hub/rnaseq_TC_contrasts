# Functions for GO Enrichment Analysis, including generating figures

gen_GO_group <- function(go_id, sig_go_df){
  # Get all the offspring GO IDs of a target ID in a list of significant
  #  GO IDs
  # Requires topGO library
  # INPUTS #
  ## go_id: target GO ID; ex: GO:0048316
  ## sig_go_df: data.frame that includes a column titled 'GO.ID' that contains
  ##  all the significant GO IDs
  # OUTPUT #
  ## list of two vectors. list[[1]] is vector containing the target GO ID and
  ##  all offspring IDs found in sig_go_df$'GO.ID'; list[[2]] is vector of same
  ##  length as list[[1]] is target GO ID (go_id) repeated - is used for making
  ##  data.frame for plotting
  ############### 
  top_id <- go_id
  kids <- GOBPOFFSPRING[[top_id]]
  kid_sig <- intersect(kids, sig_go_df$'GO.ID')
  tot_sig <- c(top_id, kid_sig)
  top_vec <- rep(top_id, times = length(tot_sig))
  go_sig_list <- list(tot_sig, top_vec)
  return(go_sig_list)
}

#######
# TEST #
### library(topGO)
### gen_go_enrich_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_general_DE_GO_enrich.txt'
### gen_go_enrich <- read.table(gen_go_enrich_file, header = T,
###   stringsAsFactors = F, sep = '\t')
### seed_dev_list <- gen_GO_group(go_term = 'GO:0048316', 
###  sig_go_df = gen_go_enrich)
########

gen_GO_barplot_df <- function(target_id_vec, sig_go_df, sig_cat_cutoffs){
  # Generate data.frame with info about GO terms to generate barplots that
  #  can be separated into different groups based on target_id_vec
  # Requires topGO package and the gen_GO_group() function from above
  # INPUTS #
  ## target_id_vec: vector of GO IDs to be processed with gen_GO_group to find
  ##  significant offspring GO IDs; the GO IDs should be high up in the
  ##  hierarchy and/or representative of a group of GO terms
  ## sig_go_df: data.frame that includes a columns titled 'GO.ID', 'Term', 
  ##  'Significant', and 'classic
  ## sig_cat_cutoffs: vector that contains the lower p-value for the different
  ##  significance categories
  # OUTPUTS #
  ## data.frame with columns to be used for making barplots in ggplot2
  ################
  tmp_list <- lapply(target_id_vec, gen_GO_group, sig_go_df = sig_go_df)
  go_id_vec <- unlist(lapply(tmp_list, function(x) x[[1]]))
  top_id_vec <- unlist(lapply(tmp_list, function(x) x[[2]]))
  tmp_df <- data.frame(go_id = go_id_vec, top = top_id_vec,
    stringsAsFactors = F)
  #
  tmp_df$term <- NA
  tmp_df$n_sig <- NA
  tmp_df$p_val <- NA
  for(i in seq(nrow(tmp_df))){
    tmp_ind <- which(sig_go_df$'GO.ID' == tmp_df$go_id[i])
    tmp_df[i, c('term', 'n_sig','p_val')] <- sig_go_df[tmp_ind, c('Term',
      'Significant', 'classic')]
  }
  #
  tmp_df$p_val[grep('< 1', tmp_df$p_val, fixed = T)] <- '1e-30'
  tmp_df$p_val <- as.numeric(tmp_df$p_val)
  tmp_df$term_go_combo <- paste(tmp_df$term, '\n', tmp_df$go_id, sep = '')
  tmp_df$top <- factor(tmp_df$top, levels = target_id_vec)
  #assign p-values to categories for plotting
  tmp_df$sig_cat <- NA
  for(psc in seq(length(sig_cat_cutoffs))){
    cat_inds <- which(tmp_df$p_val < sig_cat_cut[psc])
    tmp_df$sig_cat[cat_inds] <- psc
  }
  tmp_df$sig_cat <- factor(tmp_df$sig_cat,
    levels = sort(unique(tmp_df$sig_cat), decreasing = T))
  return(tmp_df)
}

##########
# TEST # 
### library(topGO)
### gen_go_enrich_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_general_DE_GO_enrich.txt'
### gen_go_enrich <- read.table(gen_go_enrich_file, header = T,
###   stringsAsFactors = F, sep = '\t')
### gen_target_goIDs <- c('GO:0048316', 'GO:0010876', 'GO:0006629', 
###   'GO:0098542', 'GO:0031640')
### sig_cat_cut <- c(1e-2, 1e-3, 1e-5, 1e-10)
### gen_df_long <- gen_GO_barplot_df(target_id_vec = gen_target_goIDs,
###  sig_go_df = gen_go_enrich, sig_cat_cutoffs = sig_cat_cut)
##############

# quit(save = 'no')

