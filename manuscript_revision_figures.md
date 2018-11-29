# Notes about Generating Figures Requested by Manuscript Reviewer

## Overview
* Reviewer asked for a couple figures related to the RNA-seq analysis:
  * PCA plot of the RNA-seq libraries
  * Heatmap of expression of RNA-seq DE genes

## RNA-seq PCA
### Summary
* Followed Bioconductor RNA-seq vignette to generate PCA figure
  * `http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples`
* I generated figures using two types of transformations of the read count \
data, Variance Stabilizing Transformation (VST) and rlog, but both produce \
very similar results, so I decided to use the VST data since its faster to
generate
### Workflow
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/PCA_fig.r`
### Location of Figure
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/PCA_vst.pdf`

## RNA-seq Heatmaps
### Summary
* Followed Bioconductor RNA-seq vignette to generate heatmaps of expression /
levels of the DE genes
  * `http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix`
* I generated expression heatmaps for each of the 3 DE gene lists
* I also tried generating a log2-change heatmap for the time-course /
comparisons using a Bioconductor vignette, but abandoned it. If want to \
pursue that, the codes at the end of the script after the quit() command
  * `http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments`
### Workflow
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/heatmap_figs.r`
### Location of Figures
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/GeneList1_heatmap.pdf`
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/GeneList2_heatmap.pdf`
* `/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/manuscript_figs/GeneList3_heatmap.pdf`


