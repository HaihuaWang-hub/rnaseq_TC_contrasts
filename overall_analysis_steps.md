# Steps for Camelina seed Time Course Experiment

## Data filtering
### Remove Bad Samples from Count Data
Two libraries cluster together rather than with their same-sample/same-time 
replicate libraries - data from these libraries is removed from the count data 
files
```
Rscript --vanilla /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/remove_bad_samp_counts.r
```
## Homology to Arabidopsis and GO Assignment
### Goals
* Identify top homologous Arabidopsis genes based on peptide allignment
  * At genes important when making lists of DE genes for Cs
* Assign GO terms to Camelina genes based on GO terms of top Arabidopsis
homolog
  * Will be used for GO enrichment
### Notes about Homology and GO Assignment
* On GitHub
  * https://github.com/grabowsp/rnaseq_TC_contrasts/blob/master/At_homology_and_GO_assignment.md
* On cluster
  * /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/At_homology_and_GO_assignment.md

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
### Generate Lists of DE Genes
#### Overview
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
#### Generate Gene Lists
##### Script to Make List 1: General DE
* /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/make_General_DE_List1.r
##### Script to Make List 2: Time-specific DE Genes
* /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/make_timespecific_DE_List2.r
##### Script to Make List 3: Different Response Across Time
* /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/make_different_response_DE_List3.r
#### Location of Gene Lists
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

## GO Enrichment
### Goals
* Find enriched GO Terms in each of the DE Gene Lists
* Get genes for enriched GO Terms related to lipid metabolism and seed 
development
* Generate figures to summarize/explain the GO enrichment results
### Summary
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/GO_enrichment.md`
### Workflow Scripts and Files
#### GO Enrichment Script
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/DE_GO_enrichment.r`
#### GO Enrichment Analysis and Figures Script
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/GO_result_analysis.r`
##### Custom R Functions for GO Analysis and Figures
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_tools/GO_analysis_functions.r` 

### GO Term Enrichment Lists
#### General DE GO Enrichment
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_general_DE_GO_enrich.txt`
#### DESeq2 Time Course GO Enrichment
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_TC_GO_enrich.txt`
#### splineTimeR Time-related GO Enrichment
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_TC_GO_enrich.txt`

### GO Enrichment Barplots
#### General DE GO Enrichement Barplot
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/gen_results_GO_barplot.pdf`
#### DESeq2 Time Course GO Enrichment Barplot
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/time_results_GO_barplot.pdf`
#### splineTimeR Time-related GO Enrichement Barplot
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/spline_results_GO_barplot.pdf`

## Figures for Manuscript Revision
### Summary
* Reviewer asked for a couple figures related to the RNA-seq analysis:
  * PCA plot of the RNA-seq libraries
  * Heatmap of expression of RNA-seq DE genes
### Notes
* `/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/manuscript_revision_figures.md`
### Figures
#### PCA
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/PCA_vst.pdf`
#### Heatmaps
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/GeneList1_heatmap.pdf`
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/GeneList2_heatmap.pdf`
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/GeneList3_heatmap.pdf`


