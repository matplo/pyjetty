#! /bin/bash

# This script takes a skimmed hdf5 file path as an argument, and runs a script to 
# compute nsubjettiness values for the jets/
# The main use is to give this script to a slurm script.

# Take two command line arguments -- (1) input file path, (2) output dir prefix
if [ "$1" != "" ]; then
  INPUT_FILE=$1
  #echo "Input file: $INPUT_FILE"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  JOB_ID=$2
  echo "Job ID: $JOB_ID"
else 
  echo "Wrong command line arguments"
fi

if [ "$3" != "" ]; then
  TASK_ID=$3
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

PREFIX="/rstorage/ml/egml/nsubjettiness/$JOB_ID"
SUFFIX=$(echo $INPUT_FILE | cut -d/ -f6-7)
echo $SUFFIX
OUTPUT_DIR="$PREFIX/output/$SUFFIX/$TASK_ID"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

LOG_DIR="$PREFIX/logs"
mkdir -p $LOG_DIR

# Load modules
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list

# Run python script via pipenv
cd /software/users/james/pyjetty/pyjetty/alice_analysis
pipenv run python process/user/ml/process_ppAA.py -f $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mv /rstorage/ml/egml/nsubjettiness/slurm-${JOB_ID}_${TASK_ID}.out $LOG_DIR

# Create file list
find "$PREFIX/output" -name "*.h5" >"$PREFIX/files.txt"
