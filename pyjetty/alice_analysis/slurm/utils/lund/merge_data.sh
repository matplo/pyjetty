#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=1227697

FILE_DIR="/rstorage/alice/AnalysisResults/lund/$JOB_ID"
FILES=$( find "$FILE_DIR" -name "*.root" )
echo "Number of files: $(wc -l $FILES)"

OUTPUT_DIR=/rstorage/alice/AnalysisResults/lund/$JOB_ID
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
