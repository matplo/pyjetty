#!/bin/bash

FILE_PATHS='/rstorage/u/alice/LHC18qr/147-148/files.txt'
NFILES=$(cat $FILE_PATHS | wc -l)
echo "N files to process: ${NFILES}"

STOP=$(( NFILES ))
# STOP=$(( NFILES / 10 ))
# STOP=$(( 2 ))
START=$(( 1 ))

if (( $STOP > $NFILES )) ; then
  STOP=$NFILES
fi

echo "START=$START"
echo "STOP=$STOP"

for (( JOB_N = $START; JOB_N <= $STOP; JOB_N++ ))
do
	INFILE=$(sed "${JOB_N}q;d" ${FILE_PATHS})
	OUTPUT_DIR="/storage/u/ploskon/PbPb"
	mkdir -p ${OUTPUT_DIR}
	OUTPUT_PREFIX="${OUTPUT_DIR}/${JOB_N}_rg_PbPb"
	echo ${OUTPUT_PREFIX}
	sbatch ./process_PbPb_mat.sh $INFILE ${OUTPUT_PREFIX}
done
