#!/bin/csh
set WORK_DIR=/home/t4c1/WORK/grabowsk/data/Camelina_peptide_allign
set DI_DB=arabidopsis
set QUER_FA=Cs_genes_v2.pep.fa
set OUT_FILE=cam_arab_matches.m8

cd $WORK_DIR

/home/grabowsky/bin/diamond blastp -d $DI_DB -q $QUER_FA -o $OUT_FILE --sensitive --top 10 --query-cover 50 -b1.0 -f 6 qseqid sseqid qlen slen evalue bitscore length pident

