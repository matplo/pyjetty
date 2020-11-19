#! /bin/bash
#
# Script to merge output ROOT files
SUBDIR="ang"
JOB_ID=233906
OUTPUT_DIR="/rstorage/alice/AnalysisResults/$SUBDIR/$JOB_ID"

# Merge all output files from each pt-hat bin
NBINS=20
for BIN in $(seq 1 $NBINS);
do
  FILES=$( find /rstorage/alice/AnalysisResults/$SUBDIR/$JOB_ID/449/child_*/TrainOutput/$BIN -name "*.root" )

  mkdir -p $OUTPUT_DIR/Stage1/$BIN
  # -T in below command ignores all trees
  hadd -T -f $OUTPUT_DIR/Stage1/$BIN/AnalysisResults.root $FILES

done

