#! /bin/bash
#
# Script to merge output ROOT files
JOB_ID=31173
OUTPUT_DIR="/remote_storage/hiccup6/u/alice/AnalysisResults/$JOB_ID"

# Merge all output files from each pt-hat bin
NBINS=20
for BIN in $(seq 1 $NBINS);
do

    # First, merge all output files on each hiccup
    HICCUPS="1 6 7 8 9 10 12 13"
    for HICCUP in $HICCUPS
    do
	FILE_DIR_BASE="/remote_storage/hiccup$HICCUP/u/alice/AnalysisResults/$JOB_ID"
	FILES=$( find $FILE_DIR_BASE/LHC18b8/146/child_*/TrainOutput/*/$BIN -name "*.root" )

	mkdir -p $FILE_DIR_BASE/Stage0/$BIN
	hadd -f $FILE_DIR_BASE/Stage0/$BIN/AnalysisResults_hiccup$HICCUP.root $FILES
    done

    # Then, merge each hiccup output file together into a final output file (for each bin)
    OUTPUT_DIR="/remote_storage/hiccup6/u/alice/AnalysisResults/$JOB_ID/Stage1/$BIN"
    mkdir -p $OUTPUT_DIR
    hadd -f $OUTPUT_DIR/AnalysisResults.root /remote_storage/hiccup*/u/alice/AnalysisResults/$JOB_ID/Stage0/$BIN/*.root

done
