#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=493675

OUTPUT_DIR=/rstorage/alice/AnalysisResults/james/$JOB_ID

FILE_LIST1=$OUTPUT_DIR/files_half1.txt
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal_half1.root @$FILE_LIST1

FILE_LIST2=$OUTPUT_DIR/files_half2.txt
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal_half2.root @$FILE_LIST2
