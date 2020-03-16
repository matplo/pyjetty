#! /bin/bash

# Script to merge output ROOT files
JOB_ID=40175
OUTPUT_DIR="/remote_storage/hiccup6/u/alice/AnalysisResults/$JOB_ID"

# command line arguments
if [ "$1" != "" ]; then
  MERGE_JOB_ID=$1
  echo "Merge Job ID: $MERGE_JOB_ID"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  BIN=$2
  echo "Bin: $BIN"
else
  echo "Wrong command line arguments"
fi

# Load modules
module use /software/users/james/heppy/modules
module load heppy/main_python
module use /software/users/james/pyjetty/modules
module load pyjetty/main_python
module list

# Merge all output files from each pt-hat bin
FILE_DIR_BASE=/rstorage/u/alice/AnalysisResults/$JOB_ID
FILES=$( find ${FILE_DIR_BASE}/LHC18b8/146/child_*/TrainOutput/*/${BIN} -name "*.root" )

OUT_DIR_BASE=/storage/u/alice/AnalysisResults/$JOB_ID
mkdir -p ${OUT_DIR_BASE}/Stage0/${BIN}
hadd -f -j 10 ${OUT_DIR_BASE}/Stage0/${BIN}/AnalysisResults.root $FILES

# Move stdout to appropriate folder
mv /storage/u/alice/AnalysisResults/slurm-${MERGE_JOB_ID}_${BIN}.out /storage/u/alice/AnalysisResults/${JOB_ID}/
