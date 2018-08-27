# Metadata management for Camelina seed RNAseq timecourse
## Remove Bad Samples and Add Columns Used by DE packages
Initial table is adapted from library info provided by JGI
```
samp_info_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_Timecourse_Sample_Info.tsv'
samp_info <- read.table(samp_info_file, header = T, stringsAsFactors = F)

samp_info$Treatment <- unlist(lapply( strsplit(samp_info$sampleName, 
  split = '-') , function(x) x[[1]]))
samp_info$Replicate <- unlist(lapply( strsplit(samp_info$sampleName, 
  split = '-') , function(x) x[[3]]))
samp_info$Time <- as.numeric( gsub('D', '', lapply( 
  strsplit(samp_info$sampleName, split = '-') , function(x) x[[2]]) ) )
samp_info$SampleName <- samp_info$Lib

# QC by JGI found two samples with different conditions that sample closely 
#   together - these should be removed
bad_samps <- which(samp_info$conditionNumber == 0)
samp_info_1 <- samp_info[-bad_samps, ]

metadata_file <- '/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt'
write.table(samp_info_1, file = metadata_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

```
