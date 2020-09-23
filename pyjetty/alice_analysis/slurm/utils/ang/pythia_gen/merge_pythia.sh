#! /bin/bash

# Script to merge output ROOT files
JOB_ID=189878
OUTPUT_DIR="/rstorage/alice/AnalysisResults/ang/$JOB_ID"

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
module use /home/ezra/heppy/modules
module load heppy/main_python
module use /home/ezra/pyjetty/modules
module load pyjetty/main_python
module list

# Merge all output files from each pt-hat bin
FILE_DIR_BASE=/rstorage/alice/AnalysisResults/ang/$JOB_ID
FILES=$( find ${FILE_DIR_BASE}/${BIN}/*/ -name "*.root" )

OUT_DIR_BASE=/rstorage/alice/AnalysisResults/ang/$JOB_ID
mkdir -p ${OUT_DIR_BASE}/Stage0/${BIN}
hadd -f -j 10 ${OUT_DIR_BASE}/Stage0/${BIN}/AnalysisResults.root $FILES

# Move stdout to appropriate folder
mv /rstorage/alice/AnalysisResults/ang/slurm-${MERGE_JOB_ID}_${BIN}.out /rstorage/alice/AnalysisResults/ang/${JOB_ID}/
