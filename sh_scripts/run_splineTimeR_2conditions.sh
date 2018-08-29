#!/bin/csh
set WORK_DIR=/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/sh_scripts

set SCRIPT_FILE=/home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/splineTimeR_2conditions_script.r
set READCOUNT_FILE=/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/counts_edited_v1.0.txt
set META_FILE=/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Camelina_TC_metadata_v1.0.txt
set OUT_PRE=Camelina_TC
set CONTROL_NAME=MT5

cd $WORK_DIR

Rscript --vanilla $SCRIPT_FILE $READCOUNT_FILE $META_FILE $OUT_PRE $CONTROL_NAME

