#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=1005381
OUTPUT_DIR=/rstorage/alice/AnalysisResults/ang/$JOB_ID

# Merge separate subsets, since otherwise it is too large for hadd
RUNLIST_LHC18q="000296414 000296510 000296379 000296377 000296309 000296433 000296068 000296133 000296423 000296065 000296550 000295588 000295586 000296270"
for RUN in $RUNLIST_LHC18q
do
   FILE_DIR=$OUTPUT_DIR/LHC18q/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC18r="000296894 000297446 000297544 000296899 000297479 000297442 000297415 000296934 000297590 000297380 000297123 000296694 000296903 000297218"
for RUN in $RUNLIST_LHC18r
do
   FILE_DIR=$OUTPUT_DIR/LHC18r/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

FILES=$( find $OUTPUT_DIR/LHC18*/ -name "AnalysisResultsIntermediate.root" )
hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
