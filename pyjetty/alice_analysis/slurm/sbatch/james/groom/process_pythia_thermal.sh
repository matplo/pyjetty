#! /bin/bash

# This script takes an input file path as an argument, and runs a python script to 
# process the input file and write an output ROOT file.
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

OUTPUT_DIR="/rstorage/james/groom/AnalysisResults/$JOB_ID/$TASK_ID"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
module use /software/users/james/heppy/modules
module load heppy/main_python
module use /software/users/james/pyjetty/modules
module load pyjetty/main_python
module list

# Run python script via pipenv
cd /software/users/james/pyjetty/pyjetty/alice_analysis
pipenv run python process/user/james/process_groomers.py -c config/theta_g/PbPb/james_groomers_thermal.yaml -f $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mv /rstorage/james/groom/AnalysisResults/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/james/groom/AnalysisResults/${JOB_ID}/
