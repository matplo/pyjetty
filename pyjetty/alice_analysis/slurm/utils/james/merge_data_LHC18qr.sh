#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=76667
OUTPUT_DIR=/rstorage/alice/AnalysisResults/james/$JOB_ID

# Merge separate subsets, since otherwise it is too large for hadd
RUNLIST_LHC18q="000296419 000296549 000296377 000296244 000296194 000296433 000296196 000295831 000296551 000296550 000296246 000295943 000295854 000296241"
for RUN in $RUNLIST_LHC18q
do
   FILE_DIR=$OUTPUT_DIR/LHC18qr/413-414/413/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC18r="000297315 000297414 000296935 000296850 000297544 000297481 000297442 000297441 000296690 000296941 000297222 000296934 000296852 000297379"
for RUN in $RUNLIST_LHC18r
do
   FILE_DIR=$OUTPUT_DIR/LHC18qr/413-414/414/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

FILES=$( find $OUTPUT_DIR/LHC18qr/413-414 -name "AnalysisResultsIntermediate.root" )
hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
