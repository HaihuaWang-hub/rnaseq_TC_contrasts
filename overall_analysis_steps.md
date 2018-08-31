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
  * https://github.com/grabowsp/rnaseq_TC_contrasts/blob/master/At_homology_and_GO_annotation.md
* On cluster
  * /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/At_homology_and_GO_annotation.md
## Generate Lists of DE Genes
### Overview
Generate 3 DE gene lists that fall into two categories:
1. Overall/global DE - DE across entire experiment
  * List 1: General DE
    * Overlap in DE genes from DESeq2-General, maSigPro, and 
splineTimerR-General
2. Time-Related DE
  * List 2: Time-specific DE Genes
    * Sample x Time interactions from DESeq2
  * List 3: Different Changes in Expression Through Time
    * splineTimeR spline-related results
### Generate Gene Lists
#### Script to Make List 1: General DE
* /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/make_General_DE_List1.r
#### Script to Make List 2: Time-specific DE Genes
* /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/make_timespecific_DE_List2.r
#### Script to Make List 3: Different Response Across Time
* /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/make_different_response_DE_List3.r
### Location of Gene Lists
1. List 1: General DE
  * /home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList1_DE_genes_v2.0.txt
2. List 2: Time-specific DE
  * /home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList2_DE_genes_v2.0.txt
3. List 3: Different Response Across Time DE
  * /home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList3_DE_genes_v2.0.txt
### Stats About Gene Lists
#### List 1: General DE
* 5,391 Genes; overlap between DESeq2-General, maSigPro, and
splineTimerR-General
  * 21,494 in DESeq2-General
  * 10,576 in maSigPro
  * 7,449 in splineTimeR-General
* 4,415 (81.9%) not found in any of the time-related gene sets
#### List 2: Time-Specific DE Genes
* 3,274 Genes
  * 695 (21.2%) overlap with Gene Set 1
  * 101 (3.1%) overlap with Gene Set 3
  * 2,532 (77.3%) not found in other gene sets
#### List 3: Different Response Across Time
* 433 Genes
  * 335 (77.4%) overlap with Gene Set 1
  * 101 (23.3%) overlap with Gene Set 2
  * 51 (11.8%) not found in other gene sets

