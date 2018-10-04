# Script to provide additional info about GO enrichment for Results section

# LOAD PACKAGES #
library(ggplot2)
library(topGO)
library(RColorBrewer)

# LOAD FILES/DATA #
gen_go_enrich_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_general_DE_GO_enrich.txt'
gen_go_enrich <- read.table(gen_go_enrich_file, header = T, 
  stringsAsFactors = F, sep = '\t')

time_go_enrich_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_TC_GO_enrich.txt'
time_go_enrich <- read.table(time_go_enrich_file, header = T,
  stringsAsFactors = F, sep = '\t')

spline_go_enrich_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_TC_GO_enrich.txt'
spline_go_enrich <- read.table(spline_go_enrich_file, header = T,
  stringsAsFactors = F, sep = '\t')

cs_go_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_GO_table.txt'
geneID2GO <- readMappings(file = cs_go_file)
tot_names <- names(geneID2GO)

# SET CONSTANTS #
## the colors from a palette, depending on how many sig categories there are
color_intervals <- list(6, c(5,8), c(4,7,9), c(4,6,8,9), c(4,6,7,8,9), c(4:9))

# SET VARIABLES #
## lower-cutoffs for discrete p-value categories 
sig_cat_cut <- c(1e-2, 1e-3, 1e-5, 1e-10)

## labels for the legend based on the significance cutoffs
sig_leg_labels <- rev(c('10^-2 to 10^-3', '10^-3 to 10^-5', '10^-5 to 10^-10', 
  '< 10^-10'))
## use rev() because last label should be top line in legend

# set colors for barplots
sig_col_palette <- 'Blues'
sig_colors <- brewer.pal(9, sig_col_palette)[color_intervals[[
   length(sig_cat_cut)]]]

# SET OUTPUT INFO #
fig_dir <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/'

dev_test_fig <- 'seed_dev_test_barplot.png'
gen_results_test_fig <- 'gen_results_test_barplot.png'
gen_GO_barplot_pdf <- 'gen_results_GO_barplot.pdf'

time_GO_barplot_pdf <- 'time_results_GO_barplot.pdf'
############
# topGO sandbox
### seed_dev_GO <- 'GO:0048316'
### seed_dev_kids <- GOBPOFFSPRING[[seed_dev_GO]]
### dev_kid_sig <- intersect(seed_dev_kids, gen_go_enrich$'GO.ID')
### seed_dev_tot <- c(seed_dev_GO, dev_kid_sig)

# FUNCTIONS #

# find offspring GO ids in list of signifcant GO ids
gen_GO_group <- function(go_id, sig_go_df){
  top_id <- go_id
  kids <- GOBPOFFSPRING[[top_id]]
  kid_sig <- intersect(kids, sig_go_df$'GO.ID')
  tot_sig <- c(top_id, kid_sig)
  top_vec <- rep(top_id, times = length(tot_sig))
  go_sig_list <- list(tot_sig, top_vec)
  return(go_sig_list)
}

# test with seed development
### seed_dev_list <- gen_GO_group(go_term = 'GO:0048316', 
###  sig_go_df = gen_go_enrich)
########

# function to generate data.frame with info for ggplot2 barplot
gen_GO_barplot_df <- function(target_id_vec, sig_go_df, sig_cat_cutoffs){
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

# GENERAL RESULTS BARPLOTS

## general result parent terms
### seed development: GO:0048316
### lipid localization: GO:0010876
### lipid metabolic process: GO:0006629
### defense response to other organism: GO:0098542
### cell killing: GO:0031640
gen_target_goIDs <- c('GO:0048316', 'GO:0010876', 'GO:0006629', 'GO:0098542', 
  'GO:0031640')

gen_df_long <- gen_GO_barplot_df(target_id_vec = gen_target_goIDs, 
  sig_go_df = gen_go_enrich, sig_cat_cutoffs = sig_cat_cut)

# Remove redundant GO ids
## some GO IDs contain the exact same sets of genes and have VERY similar
##  terms, so removing them a) makes the figure cleaner and b) is a more
##  representative display of the results because it does not over-represent
##  certain data
gen_ids_remove <- c('GO:0046836', 'GO:0033559', 'GO:0033383', 'GO:0033385')
gen_remove_inds <- c()
for(gir in gen_ids_remove){
  tmp_ind <- which(gen_df_long$go_id == gir)
  gen_remove_inds <- c(gen_remove_inds, tmp_ind)
}

gen_df <- gen_df_long[-gen_remove_inds, ]

# change facet text formatting so fits in "strip"
gen_facet_labs <- sapply(gen_target_goIDs, Term)
gen_facet_labs[grep('defense response', 
  gen_facet_labs, fixed = T)] <- 'defense response\nto other organism'
gen_facet_labs[grep('killing', gen_facet_labs)] <- ' '
names(gen_facet_labs) <- gen_target_goIDs

genresults_barplot <- ggplot(gen_df, aes(x = reorder(term_go_combo, -n_sig), 
  y = n_sig, fill = sig_cat)) +
  # the 'reorder' command orders the bars by n_sig
  facet_grid(~top, scales = 'free_x', space = 'free_x', 
    labeller = as_labeller(gen_facet_labs)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(trans = 'log2') +
  ylab('DE genes in GO Category (log2 scale)') +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(name = 'p-value', breaks = c(4,3,2,1), 
    labels = sig_leg_labels, values = sig_colors) 

#png(filename = paste(fig_dir, gen_results_test_fig, sep = ''))
#genresults_barplot
#dev.off()

pdf(file = paste(fig_dir, gen_GO_barplot_pdf, sep = ''), width = 9, height = 9)
genresults_barplot
dev.off()

# TIME-SPECIFIC RESULTS BARPLOTS

## time-specific result parent terms
### seed development: GO:0048316
### lipid localization: GO:0010876
### lipid metabolic process: GO:0006629
### photosynthesis: GO:0015979
### defense response to other organism: GO:0098542
time_target_goIDs <- c('GO:0048316', 'GO:0010876', 'GO:0006629', 'GO:0015979', 
  'GO:0098542')

time_df_long <- gen_GO_barplot_df(target_id_vec = time_target_goIDs,
  sig_go_df = time_go_enrich, sig_cat_cutoffs = sig_cat_cut)

# Remove redundant GO ids
time_ids_remove <- c('GO:0033559', 'GO:0016131', 'GO:0016128', 'GO:0016129', 'GO:0016103', 'GO:0043155', 'GO:1905156', 'GO:0019685')
time_remove_inds <- c()
for(ir in time_ids_remove){
  tmp_ind <- which(time_df_long$go_id == ir)
  time_remove_inds <- c(time_remove_inds, tmp_ind)
}

time_df <- time_df_long[-time_remove_inds, ]

# change facet text formatting so fits in "strip"
time_facet_labs <- sapply(time_target_goIDs, Term)
time_facet_labs[grep('defense response',
  time_facet_labs, fixed = T)] <- 'defense\nresponse to\nother organism'
time_facet_labs[grep('localization', time_facet_labs)] <- 'lipid\nlocalization'
names(time_facet_labs) <- time_target_goIDs

time_results_barplot <- ggplot(time_df, aes(x = reorder(term_go_combo, -n_sig),
  y = n_sig, fill = sig_cat)) +
  # the 'reorder' command orders the bars by n_sig
  facet_grid(~top, scales = 'free_x', space = 'free_x',
    labeller = as_labeller(time_facet_labs)) +
  geom_col(position = "dodge") +
  scale_y_continuous(trans = 'log2') +
  ylab('DE genes in GO Category (log2 scale)') +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(name = 'p-value', breaks = c(4,3,2,1),
    labels = sig_leg_labels, values = sig_colors)

#png(filename = paste(fig_dir, gen_results_test_fig, sep = ''))
#genresults_barplot
#dev.off()

pdf(file = paste(fig_dir, time_GO_barplot_pdf, sep = ''), 
  width = 15, height = 9)
time_results_barplot
dev.off()


# SPLINE RESULTS BOXPLOT
# THIS IS COPIED FROM ABOVE - NEED TO ADJUST
## time-specific result parent terms
### seed development: GO:0048316
### lipid localization: GO:0010876
### lipid metabolic process: GO:0006629
### photosynthesis: GO:0015979
### defense response to other organism: GO:0098542
time_target_goIDs <- c('GO:0048316', 'GO:0010876', 'GO:0006629', 'GO:0015979',
  'GO:0098542')



########

quit(save = 'no')

