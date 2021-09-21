#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=611758

FILE_DIR="/rstorage/alice/AnalysisResults/ang/$JOB_ID"
FILES=$( find "$FILE_DIR" -name "*.root" )

OUTPUT_DIR=/rstorage/alice/AnalysisResults/ang/$JOB_ID
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
