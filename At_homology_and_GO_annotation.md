# Steps for assigning homologous Arabidopsis thaliana genes and GO terms to
Camelina genes
## Goals
1) Align Camelina peptides to Arabidopsis pepetides
2) Assign top Arabidopsis matches from alignments to Camelina genes
3) Assign GO terms to Camelina genes using GO terms assigned to homologous
Arabidopsis genes
## Align Camelina Pepetide Sequences to Arabidopsis
### Resources
* I use the program `diamond` to quickly use BLASTP for many genes
* The file `Araport11_genes.201606.pep.fasta` was downloaded from TAIR and used
for the Arabidopsis peptide sequences
* The file `Cs_genes_v2.pep.fa` was downloaded from http://www.camelinadb.ca/downloads.html
### Workflow
#### Generate Reference Database
Required by diamond
```
cd /home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign
diamond makedb --in Araport11_genes.201606.pep.fasta -d arabidopsis
```
#### Test diamond on Small Set of Genes
```
cd /home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign
diamond blastp -d arabidopsis -q test_Cs.pep.fa -o test_matches.m8 --sensitive --top 10 --query-cover 50 -b1.0 -f 6 qseqid sseqid qlen slen evalue bitscore length pident
# next step is doing the alignment
diamond blastp -d arabidopsis -q READS.FA -o cam_arab_matches.m8 --sensitive --top 10 --query-cover 50 -b1.0 -f 6 qseqid sseqid qlen slen evalue bitscore length pident
```
#### Full Query of Camelina vs Arabidopsis peptides
##### Shell script
`/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/sh_scripts/run_diamond_Cam_peps.sh`
##### Submit Full Queery Job
```
qsub -cwd -N cam_arab_blastp -l h_vmem=7G -q all.q run_diamond_Cam_peps.sh
```
## Assign A.thaliana homologs and GO terms to Camelina Genes
### Resources
* Arabidopsis GO terms in file `ATH_GO_GOSLIM.txt` downloaded from TAIR
* Arabidopsis gene alias file `gene_aliases_20130831.txt` downloaded from TAIR
### Workflow
#### Assign homologs and GO terms
```
Rscript --vanilla /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/assign_At_homologs_and_GO.r
```
#### Assign Arabidopsis common name to A.thaliana homologs
```
Rscript --vanilla /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/assign_At_common_name.r
```

