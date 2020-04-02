#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=44195

# First, merge all output files on each hiccup
FILE_DIR="/rstorage/alice/AnalysisResults/ang/$JOB_ID"
FILES=$( find "$FILE_DIR" -name "*.root" )

mkdir -p $FILE_DIR/Stage1
hadd -f $FILE_DIR/Stage1/AnalysisResults_hiccup$HICCUP.root $FILES

# Then, merge each hiccup output file together into a final output file
OUTPUT_DIR="/rstorage/alice/AnalysisResults/ang/$JOB_ID"
hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root /rstorage/alice/AnalysisResults/ang/$JOB_ID/Stage1/*.root
