#! /bin/bash
#
# Script to merge output ROOT files from all pt-hat bins together, in stages
JOB_ID=61683
FILE_DIR=/rstorage/alice/AnalysisResults/ang/$JOB_ID
OUTPUT_DIR=/rstorage/alice/AnalysisResults/ang/$JOB_ID

# Merge all output files from each pt-hat bin
mkdir -p $OUTPUT_DIR
hadd -f -j 10 $OUTPUT_DIR/AnalysisResultsFinal.root $FILE_DIR/Stage0/*/*.root
