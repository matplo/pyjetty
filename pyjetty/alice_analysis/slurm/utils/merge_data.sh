#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=40084

FILE_DIR="/rstorage/u/alice/AnalysisResults/$JOB_ID"
FILES=$( find "$FILE_DIR" -name "*.root" )

OUTPUT_DIR=/remote_storage/hiccup6/u/alice/AnalysisResults/$JOB_ID
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
