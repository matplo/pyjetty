#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=19626

# First, merge all output files on each hiccup
HICCUPS="1 6 7 8 9 10 12 13"
for HICCUP in $HICCUPS
do
  FILE_DIR="/remote_storage/hiccup$HICCUP/u/alice/AnalysisResults/$JOB_ID"
  FILES=$( find "$FILE_DIR" -name "*.root" )

  mkdir -p $FILE_DIR/Stage1
  hadd -f $FILE_DIR/Stage1/AnalysisResults_hiccup$HICCUP.root $FILES

done

# Then, merge each hiccup output file together into a final output file
OUTPUT_DIR="/remote_storage/hiccup6/u/alice/AnalysisResults/$JOB_ID"
hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root /remote_storage/hiccup*/u/alice/AnalysisResults/$JOB_ID/Stage1/*.root
