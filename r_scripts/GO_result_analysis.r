# Script to provide additional info about GO enrichment for Results section

# LOAD PACKAGES #
library(ggplot2)
library(topGO)
library(RColorBrewer)

## load functions
go_analysis_funcs_file <- '/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_tools/GO_analysis_functions.r'
source(go_analysis_funcs_file)

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
### use rev() because last label should be top line in legend

# set colors for barplots
sig_col_palette <- 'Blues'
sig_colors <- brewer.pal(9, sig_col_palette)[color_intervals[[
   length(sig_cat_cut)]]]

# aprox. width, in inches of bars in plots
indiv_bar_width <- 0.275
# aprox. width of legend and margins
bar_extra_width <- 2.25
# all the plots will have the same height
bar_plot_height <- 9

# SET OUTPUT INFO #
fig_dir <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/'

dev_test_fig <- 'seed_dev_test_barplot.png'
gen_results_test_fig <- 'gen_results_test_barplot.png'
gen_GO_barplot_pdf <- 'gen_results_GO_barplot.pdf'

time_GO_barplot_pdf <- 'time_results_GO_barplot.pdf'

spline_GO_barplot_pdf <- 'spline_results_GO_barplot.pdf'
############
# topGO sandbox
### seed_dev_GO <- 'GO:0048316'
### seed_dev_kids <- GOBPOFFSPRING[[seed_dev_GO]]
### dev_kid_sig <- intersect(seed_dev_kids, gen_go_enrich$'GO.ID')
### seed_dev_tot <- c(seed_dev_GO, dev_kid_sig)

# OVERLAP IN GO CATEGORIES #
length(intersect(gen_go_enrich$'GO.ID', intersect(time_go_enrich$'GO.ID', 
  spline_go_enrich$'GO.ID')))
# [1] 25
# 25 GO terms are enriched in all three gene sets
#sapply(intersect(gen_go_enrich$'GO.ID', intersect(time_go_enrich$'GO.ID', 
#  spline_go_enrich$'GO.ID')), Term)

length(setdiff(gen_go_enrich$'GO.ID', union(time_go_enrich$'GO.ID', 
  spline_go_enrich$'GO.ID')))
# [1] 84
# 84 GO terms unique to List 1

length(setdiff(time_go_enrich$'GO.ID', union(gen_go_enrich$'GO.ID',      
  spline_go_enrich$'GO.ID')))
# [1] 328
# 328 GO terms unique to List 2

length(setdiff(spline_go_enrich$'GO.ID', union(gen_go_enrich$'GO.ID', 
  time_go_enrich$'GO.ID')))
# [1] 34
# 34 GO terms unique to List 3

length(intersect(time_go_enrich$'GO.ID', spline_go_enrich$'GO.ID'))
# [1] 39
# 39 GO terms shared between List 2 and List 3

length(intersect(gen_go_enrich$'GO.ID', time_go_enrich$'GO.ID'))
# [1] 77
# 77 GO terms shared between List 1 and List 2

length(intersect(gen_go_enrich$'GO.ID', spline_go_enrich$'GO.ID'))
# [1] 47
# 47 GO terms shared between List 1 and List 3

lip_met_GO <- 'GO:0006629'
lip_met_kids <- GOBPOFFSPRING[[lip_met_GO]]

length(intersect(lip_met_kids, gen_go_enrich$'GO.ID'))
# [1] 9
# 9 enriched lipid metabolism offspring terms in List 1

length(intersect(lip_met_kids, time_go_enrich$'GO.ID'))
# [1] 29
# 29 enriched lipid metabolism offspring terms in List 2

length(intersect(lip_met_kids, spline_go_enrich$'GO.ID'))
# [1] 21
# 21 enriched lipid metabolism offspring terms in List 3

###############
# GENERAL RESULTS BARPLOTS

## general result parent terms
### seed development: GO:0048316
### lipid localization: GO:0010876
### lipid metabolic process: GO:0006629
### phenylpropanoid metabolic process: GO:0009698
### defense response to other organism: GO:0098542
### killing of cells of other organism: GO:0031640
gen_target_goIDs <- c('GO:0048316', 'GO:0010876', 'GO:0006629', 'GO:0009698', 
  'GO:0098542', 'GO:0031640')

gen_df_long <- gen_GO_barplot_df(target_id_vec = gen_target_goIDs, 
  sig_go_df = gen_go_enrich, sig_cat_cutoffs = sig_cat_cut)
# gen_GO_barplot_df() is function imported from [GO_analysis_functions.r] 

# Remove redundant GO ids
## some GO IDs contain the exact same sets of genes and have VERY similar
##  terms, so removing them a) makes the figure cleaner and b) is a more
##  representative display of the results because it does not over-represent
##  certain data
gen_ids_remove <- c('GO:0046836', 'GO:0033559', 'GO:0033383', 'GO:0033385', 
  'GO:0009808', 'GO:0009806', 'GO:0009803')
gen_remove_inds <- c()
for(gir in gen_ids_remove){
  tmp_ind <- which(gen_df_long$go_id == gir)
  gen_remove_inds <- c(gen_remove_inds, tmp_ind)
}

gen_df <- gen_df_long[-gen_remove_inds, ]

# change facet text formatting so fits in "strip"
gen_facet_labs <- sapply(gen_target_goIDs, Term)
gen_facet_labs[grep('seed',
  gen_facet_labs)] <- 'seed\ndevelopment'
gen_facet_labs[grep('defense response', 
  gen_facet_labs, fixed = T)] <- 'defense\nresponse\nto other\norganism'
gen_facet_labs[grep('phenylpropanoid', 
  gen_facet_labs)] <- 'phenylpropanoid\nmetabolic process'
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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90)) +
  scale_fill_manual(name = 'p-value', breaks = c(4,3,2,1), 
    labels = sig_leg_labels, values = sig_colors) 

#png(filename = paste(fig_dir, gen_results_test_fig, sep = ''))
#genresults_barplot
#dev.off()

gen_plot_width <- nrow(gen_df)*indiv_bar_width + bar_extra_width

pdf(file = paste(fig_dir, gen_GO_barplot_pdf, sep = ''), 
  width = gen_plot_width, height = bar_plot_height)
genresults_barplot
dev.off()

#######################
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
time_ids_remove <- c('GO:0033559', 'GO:0016131', 'GO:0016128', 'GO:0016129', 
  'GO:0016103', 'GO:0043155', 'GO:1905156', 'GO:0019685')
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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90)) +
  scale_fill_manual(name = 'p-value', breaks = c(4,3,2,1),
    labels = sig_leg_labels, values = sig_colors)

time_plot_width <- nrow(time_df)*indiv_bar_width + bar_extra_width

pdf(file = paste(fig_dir, time_GO_barplot_pdf, sep = ''), 
  width = time_plot_width, height = bar_plot_height)
time_results_barplot
dev.off()

#######################
# SPLINE RESULTS BOXPLOT

## spline result parent terms
### lipid metabolic process: GO:0006629
### phenylpropanoid metaboloic process: GO:0009698 (includes suberin)
### cutin biosynthetic process: GO:0010143
### killing of cells of other orgamism: GO:0031640
spline_target_goIDs <- c('GO:0006629', 'GO:0009698', 'GO:0010143', 
  'GO:0031640')

spline_df_long <- gen_GO_barplot_df(target_id_vec = spline_target_goIDs,
  sig_go_df = spline_go_enrich, sig_cat_cutoffs = sig_cat_cut)

# Remove redundant GO ids
spline_ids_remove <- c('GO:0006644', 'GO:0006720', 'GO:0006720', 'GO:0016116',
  'GO:0016108', 'GO:0045338', 'GO:0033383', 'GO:0033385', 'GO:0000038',
  'GO:0009803')
spline_remove_inds <- c()
for(ir in spline_ids_remove){
  tmp_ind <- which(spline_df_long$go_id == ir)
  spline_remove_inds <- c(spline_remove_inds, tmp_ind)
}

spline_df <- spline_df_long[-spline_remove_inds, ]

# change facet text formatting so fits in "strip"
spline_facet_labs <- sapply(spline_target_goIDs, Term)
spline_facet_labs[grep('phenylpropanoid',
  spline_facet_labs, fixed = T)] <- 'phenylpropanoid\nmetabolic\nprocess'
spline_facet_labs[grep('cutin', spline_facet_labs)] <- ''
spline_facet_labs[grep('killing', spline_facet_labs)] <- ''
names(spline_facet_labs) <- spline_target_goIDs

spline_results_barplot <- ggplot(spline_df, 
  aes(x = reorder(term_go_combo, -n_sig), y = n_sig, fill = sig_cat)) +
  # the 'reorder' command orders the bars by n_sig
  facet_grid(~top, scales = 'free_x', space = 'free_x',
    labeller = as_labeller(spline_facet_labs)) +
  geom_col(position = "dodge") +
  scale_y_continuous(trans = 'log2') +
  ylab('DE genes in GO Category (log2 scale)') +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 90)) +
  scale_fill_manual(name = 'p-value', breaks = c(4,3,2,1),
    labels = sig_leg_labels, values = sig_colors)

spline_plot_width <- nrow(spline_df)*indiv_bar_width + bar_extra_width

pdf(file = paste(fig_dir, spline_GO_barplot_pdf, sep = ''), 
  width = spline_plot_width, height = bar_plot_height)
spline_results_barplot
dev.off()

########

quit(save = 'no')

