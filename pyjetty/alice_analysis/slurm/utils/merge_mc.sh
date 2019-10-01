#! /bin/bash
#
# Script to merge output ROOT files
JOB_ID=7467
OUTPUT_DIR="/remote_storage/hiccup6/u/alice/AnalysisResults/$JOB_ID"

# Merge all output files from each pt-hat bin
NBINS=20
for BIN in $(seq 1 $NBINS);
do
  FILES=$( find /rstorage/u/alice/AnalysisResults/$JOB_ID/LHC18b8/146/child_*/TrainOutput/*/$BIN -name "*.root" )

  mkdir -p $OUTPUT_DIR/Stage1/$BIN
  hadd -f $OUTPUT_DIR/Stage1/$BIN/AnalysisResults.root $FILES

done
