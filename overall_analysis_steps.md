# Steps for Camelina seed Time Course Experiment
## Data filtering
### Remove Bad Samples from Count Data
Two libraries cluster together rather than with their same-sample/same-time 
replicate libraries - data from these libraries is removed from the count data 
files
```
Rscript --vanilla /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/remove_bad_samp_counts.r
```
## DE Analysis
### Goals
* Identify Differentially Expressed (DE) genes using 3 R packages: DESeq2, 
maSigPro, and splineTimeR
  * Two approaches: 1) Overall DE across experiment; 2) Differential response
across time (probably the most important results given the experimental design)
* Generate lists of top DE genes that combine results from the 3 packages and
are based on the models used to identify DE genes
### Notes about DE analysis
* On GitHub
  * https://github.com/grabowsp/rnaseq_TC_contrasts/blob/master/DE_analysis_notes.md
* On cluster
  * /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/DE_analysis_notes.md
## Homology to Arabidopsis and GO Assignment
### Goals
* Identify top homologous Arabidopsis genes based on peptide allignment
* Assign GO terms to Camelina genes based on GO terms of top Arabidopsis
homolog
### Notes about Homology and GO Assignment
* On GitHub
* On cluster
  * /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/At_homology_and_GO_annotation.md
