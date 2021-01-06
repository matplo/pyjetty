#!/bin/bash

if [ "$1" != "" ]; then
  INDEX=$1
  echo "Index: $INDEX"
else
  echo "Wrong command line arguments"
fi

outputDir="/mnt/rstorage/alice/data/LHC18qr/550"

year="2018"
period="LHC18q"
trainName="HF_TreeCreator"
trainPWG="pass3/PWGHF"
trainNumber="550_20201223-0456_child_1"
filename="AnalysisResults.root"

#LHC18q pass3
RUNLIST=(000296414 000296510 000296379 000296377 000296309 000296433 000296068 000296133 000296423 000296065 000296550 000295588 000295586 000296270)
#LHC18r pass3
#RUNLIST=(000296894 000297446 000297544 000296899 000297479 000297442 000297415 000296934 000297590 000297380 000297123 000296694 000296903 000297218)

#LHC18q pass1
#RUNLIST=(000296419 000296549 000296377 000296244 000296194 000296433 000296196 000295831 000296551 000296550 000296246 000295943 000295854 000296241)
#LHC18r pass1
#RUNLIST=(000297315 000297414 000296935 000296850 000297544 000297481 000297442 000297441 000296690 000296941 000297222 000296934 000296852 000297379)
RUN=${RUNLIST[$INDEX]}

TRAIN_OUTPUT_DIR=/alice/data/${year}/${period}/$RUN/${trainPWG}/${trainName}/${trainNumber}

FILELIST=$(alien_ls $TRAIN_OUTPUT_DIR)
echo $FILELIST

MAX=$(alien_ls $TRAIN_OUTPUT_DIR | wc -l)
SUBDIRS="$(seq -w 1 $MAX)"
echo "Download $MAX files at $TRAIN_OUTPUT_DIR"

exit 1

# Create outputDir and cd into it
if [ ! -d $outputDir ]; then
  mkdir -p $outputDir
fi
cd $outputDir
echo "output dir: $outputDir"

# Remove any empty directories
if [ -d "$RUN" ]; then
  find $RUN -empty -type d -delete
fi

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
