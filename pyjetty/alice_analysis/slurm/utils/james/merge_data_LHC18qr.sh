#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=63515
OUTPUT_DIR=/rstorage/alice/AnalysisResults/james/$JOB_ID

FILE_DIR="/rstorage/alice/AnalysisResults/james/$JOB_ID/LHC18qr/147-148/147"
FILES=$( find "$FILE_DIR" -name "*.root" )
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal_147.root $FILES

FILE_DIR="/rstorage/alice/AnalysisResults/james/$JOB_ID/LHC18qr/147-148/148"
FILES=$( find "$FILE_DIR" -name "*.root" )
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal_148.root $FILES

hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root $OUTPUT_DIR/AnalysisResultsFinal_147.root $OUTPUT_DIR/AnalysisResultsFinal_148.root
