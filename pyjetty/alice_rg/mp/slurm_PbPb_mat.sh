#!/bin/bash

FILE_PATHS='/rstorage/u/alice/LHC18qr/147-148/files.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# STOP=$(( NFILES / 10 ))
STOP=$(( 2 ))
START=$(( 0 ))

if (( $STOP > $NFILES ))
then
  STOP=$NFILES
fi

echo "START=$START"
echo "STOP=$STOP"

for (( JOB_N = $START; JOB_N <= $STOP; JOB_N++ ))
do
  FILE=$(sed -n "$JOB_N"p $FILE_PATHS)
  OUTPUT_PREFIX="/storage/u/ploskon/PbPb/${JOB_N}"
  sbatch ./process_PbPb_mat.sh $FILE $OUTPUT_PREFIX
done
