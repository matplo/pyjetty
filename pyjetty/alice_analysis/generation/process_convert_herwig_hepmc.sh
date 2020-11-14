#! /bin/bash

# This script takes an input HepMC file path as an argument, and runs a python script to
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

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

# Define output path from relevant sub-path of input file
# Note: suffix depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f4-7)
echo $OUTPUT_SUFFIX
OUTPUT_DIR="/rstorage/generators/herwig_alice/$JOB_ID/$OUTPUT_SUFFIX"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list

# Run main script
python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR -g herwig --no-progress-bar

# Move stdout to appropriate folder
mv /rstorage/generators/herwig_alice/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/generators/herwig_alice/${JOB_ID}/


