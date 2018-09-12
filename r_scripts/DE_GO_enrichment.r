# Script for GO enrichment of DE genes

# LOAD FILES #
cs_go_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign/Cam_GO_table.txt'

count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')

## splineTimeR results
spline_tc_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_noIntercept_genes.txt'
spline_tc <- read.table(spline_tc_file, header = T, stringsAsFactors = F, 
  sep = '\t')

## DESeq2 time-related results
deseq_tc_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_TC_genes.txt'
deseq_tc <- read.table(deseq_tc_file, header = T, stringsAsFactors = F, 
  sep = '\t')

## Combo results
combo_gen_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_Combo_general_genes.txt'
combo_gen <- read.table(combo_gen_file, header = T, stringsAsFactors = F, 
  sep = '\t')

# SET VARIABLES #

# SET CONSTANTS #
min_floor_count_cut <- 5

# SET OUTPUT INFO #
spline_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_TC_GO_enrich.txt'

deseq_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_TC_GO_enrich.txt'

gen_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_general_DE_GO_enrich.txt'

# LOAD PACKAGES #
library('topGO', lib.loc = '/home/grabowsky/tools/r_packages')
library('edgeR', lib.loc = '/home/grabowsky/tools/r_packages')

######################
# Generate topGO object used for GO Enrichment
geneID2GO <- readMappings(file = cs_go_file)
tot_names <- names(geneID2GO)

# Filter overall gene list - only use genes include in DE analysis
counts_1 <- counts[, -1]
rownames(counts_1) <- counts[,1]

norm_factors <- calcNormFactors(counts_1, method = 'TMM')
counts_2 <- mapply(`/`, counts_1, norm_factors)
rownames(counts_2) <- rownames(counts_1)

# Remove genes with low counts across samples
max_count <- apply(counts_2, 1, max)
low_rows <- which(max_count <= min_floor_count_cut)
counts_3 <- counts_2[-low_rows, ]

expr_genes <- rownames(counts_3)
expr_inds <- which(tot_names %in% expr_genes)

geneID2GO_2 <- geneID2GO[expr_inds]
tot_names_2 <- names(geneID2GO_2)

# Function for generating GO enrichment tables from DE tables
gen_GO_enrich_df <- function(de_df, id2go_obj, sig_val = 0.01){
  go_obj_genes <- names(id2go_obj)
  test_geneList <- factor(as.integer(go_obj_genes %in% de_df$gene))
  names(test_geneList) <- go_obj_genes
  test_GOdata <- new('topGOdata', ontology = 'BP', allGenes = test_geneList, 
    annot = annFUN.gene2GO, gene2GO = id2go_obj, nodeSize = 5) 
  test.stat <- new('classicCount', testStatistic = GOFisherTest, 
    name = 'Fisher test')
  test_resultFisher <- getSigGroups(test_GOdata, test.stat)
  n_sig_terms <- sum(score(test_resultFisher) < sig_val)
  test_enrich_df <- GenTable(test_GOdata, classic = test_resultFisher, 
    orderBy = classic, ranksOf = 'classic', topNodes = n_sig_terms)
  return(test_enrich_df)
}

# Run enrichment and generate output tables
spline_tc_df <- gen_GO_enrich_df(de_df = spline_tc, id2go_obj = geneID2GO_2)
write.table(spline_tc_df, file = spline_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

deseq_tc_df <- gen_GO_enrich_df(de_df = deseq_tc, id2go_obj = geneID2GO_2)
write.table(deseq_tc_df, file = deseq_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

gen_DE_df <- gen_GO_enrich_df(de_df = combo_gen, id2go_obj = geneID2GO_2)
write.table(gen_DE_df, file = gen_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

quit(save = 'no')

