# Analysis of DE genes for manuscript

# LOAD PACKAGES #

# LOAD DATA #
list1_DE_genes_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList1_DE_genes_v2.0.txt'
list1_DE_genes <- read.table(list1_DE_genes_file, header = T, 
  stringsAsFactors = F, sep = '\t')

list2_DE_genes_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList2_DE_genes_v2.0.txt'
list2_DE_genes <- read.table(list2_DE_genes_file, header = T, 
  stringsAsFactors = F, sep = '\t')

list3_DE_genes_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList3_DE_genes_v2.0.txt'
list3_DE_genes <- read.table(list3_DE_genes_file, header = T, 
  stringsAsFactors = F, sep = '\t')

# SET CONSTANTS #


# SET VARIABLES #


# SET OUTPUT INFO #


#############

length(intersect(list1_DE_genes$gene, intersect(list2_DE_genes$gene, 
  list3_DE_genes$gene)))
# [1] 54
# 54 genes found in all 3 lists

length(setdiff(list1_DE_genes$gene, union(list2_DE_genes$gene, 
  list3_DE_genes$gene)))
# [1] 4415
# 4415 genes (81.9%) found ONLY in List 1

length(setdiff(list2_DE_genes$gene, union(list1_DE_genes$gene, 
  list3_DE_genes$gene)))
# [1] 2532
# 2532 genes (77.3%) found ONLY in List 2

length(setdiff(list3_DE_genes$gene, union(list1_DE_genes$gene,
  list2_DE_genes$gene)))
# [1] 51
# 51 genes (11.8%) found ONLY in List 3

length(intersect(list1_DE_genes$gene, list2_DE_genes$gene))
# [1] 695
# 695 genes overlap between List 1 and List 2

length(intersect(list1_DE_genes$gene, list3_DE_genes$gene))
# [1] 335
# 335 genes overlap between List 1 and List 3

length(intersect(list2_DE_genes$gene, list3_DE_genes$gene))
# [1] 101
# 101 genes overlap between List 2 and List 3

length(setdiff(intersect(list2_DE_genes$gene, list3_DE_genes$gene), 
  list1_DE_genes$gene))
# [1] 47
# 47 genes overlap between List 2 and List 3 and are NOT part of List 1



