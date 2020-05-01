#!/bin/bash

if [ "$1" != "" ]; then
  INDEX=$1
  echo "Index: $INDEX"
else
  echo "Wrong command line arguments"
fi

outputDir="/mnt/rstorage/alice/data/LHC18qr/413-414/414"

year="2018"
period="LHC18r"
trainName="HF_TreeCreator"
trainPWG="pass1/PWGHF"
trainNumber="414_20200429-0933"
filename="AnalysisResults.root"

#LHC18q
#RUNLIST=(000296419 000296549 000296377 000296244 000296194 000296433 000296196 000295831 000296551 000296550 000296246 000295943 000295854 000296241)
#LHC18r
RUNLIST=(000297315 000297414 000296935 000296850 000297544 000297481 000297442 000297441 000296690 000296941 000297222 000296934 000296852 000297379)
RUN=${RUNLIST[$INDEX]}

TRAIN_OUTPUT_DIR=/alice/data/${year}/${period}/$RUN/${trainPWG}/${trainName}/${trainNumber}

MAX=$(alien_ls $TRAIN_OUTPUT_DIR | wc -l)
SUBDIRS="$(seq -w 1 $MAX)"
echo "Download $MAX files at $TRAIN_OUTPUT_DIR"

# Create outputDir and cd into it
if [ ! -d $outputDir ]; then
  mkdir -p $outputDir
fi
cd $outputDir
echo "output dir: $outputDir"

# Remove any empty directories
find $RUN -empty -type d -delete

# Copy the files (skip any directory that already exists)
for SUBDIR in $SUBDIRS
do
  if [ ! -d "$RUN/$SUBDIR" ]; then
    mkdir -p "$RUN/$SUBDIR"
  else
    continue
  fi        
  alien_cp alien://$TRAIN_OUTPUT_DIR/$SUBDIR/AnalysisResults.root $RUN/$SUBDIR
done
