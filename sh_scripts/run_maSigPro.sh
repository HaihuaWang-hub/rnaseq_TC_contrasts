#!/bin/csh
set WORK_DIR=/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC

set SCRIPT_FILE=/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/maSigPro_TC_script.r
set READCOUNT_FILE=/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt
set META_FILE=/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt
set OUT_PRE=Camelina_TC_maSigPro
set CONTROL_NAME=MT5

cd $WORK_DIR

Rscript --vanilla $SCRIPT_FILE $READCOUNT_FILE $META_FILE $OUT_PRE $CONTROL_NAME

