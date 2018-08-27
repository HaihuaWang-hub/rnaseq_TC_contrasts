# Script to remove the bad samples from the read-count data for the Camelina
#  seed RNA-seq experiment

# LOAD FILES #
samp_info_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_Timecourse_Sample_Info.tsv'
samp_info <- read.table(samp_info_file, header = T, stringsAsFactors = F)

count_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #
count_out_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt'

# SET CONSTANTS #


# LOAD PACKAGES #


####################
meta_bad_samp_inds <- which(samp_info$conditionNumber == 0)

data_bad_inds <- c()
for(metaInd in meta_bad_samp_inds){
  tmp_ind <- which(colnames(counts) == samp_info$Lib[metaInd])
  data_bad_inds <- c(data_bad_inds, tmp_ind)
}

write.table(counts_1, file = count_out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

quit(save = 'no')
