#!/bin/bash

FILES=$(find $PWD/split100 -name "PbPbfiles_??.txt")
for fn in ${FILES}
do
	BNAME=$(basename ${fn})
	OUTPUT_DIR="/storage/u/ploskon/PbPb"
	mkdir -p ${OUTPUT_DIR}
	OUTPUT_PREFIX="${OUTPUT_DIR}/${BNAME}_rg_PbPb"
	echo ${fn} ${OUTPUT_PREFIX}
	sbatch ./process_PbPb_mat.sh ${fn} ${OUTPUT_PREFIX}
done
