# GO Enrichment Analysis of DE Genes
## Goals
1) Generate GO enrichment lists for the different DE Gene lists
2) Generate figures explaining GO enrichment results

## GO Enrichment
### Workflow
#### Generate GO Enrichment lists
```
Rscript --vanilla /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/DE_GO_enrichment.r
```
### Location of Lists
#### General DE GO Enrichment
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_general_DE_GO_enrich.txt`
#### DESeq2 Time Course GO Enrichment
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_DESeq2_TC_GO_enrich.txt`
#### splineTimeR Time-related GO Enrichment
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_splineTimeR_TC_GO_enrich.txt`

## GO Enrichment Figures and Analysis
### Overview
* Generate bargraphs showing number of DE genes in certain enriched GO \
categories
  * Grouped by higher GO terms
  * Colored by p-value
### Workflow
```
/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/GO_result_analysis.r
```
### Location of Figures
#### General DE GO Enrichement Barplot
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/gen_results_GO_barplot.pdf`
#### DESeq2 Time Course GO Enrichment Barplot
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/time_results_GO_barplot.pdf`
#### splineTimeR Time-related GO Enrichement Barplot
`/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/figs/spline_results_GO_barplot.pdf`



