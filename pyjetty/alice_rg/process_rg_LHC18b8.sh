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
  OUTPUT_PREFIX=$2
  echo "Output dir prefix: $OUTPUT_PREFIX"
else 
  echo "Wrong command line arguments"
fi

# Define output path from relevant sub-path of input file
# Note: depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f5-10)
#echo $OUTPUT_SUFFIX
OUTPUT_DIR="/storage/u/alice/$OUTPUT_PREFIX/$OUTPUT_SUFFIX"
#echo "Output dir: $OUTPUT_DIR"

# Load modules
module use /software/users/james/heppy/modules
module load heppy/main_python3
module use /software/users/james/pyjetty/modules
module load pyjetty/main_python3
module list

# Run python script via pipenv
cd /software/users/james/pyjetty/pyjetty/alice_rg
pipenv run python process_rg_data.py -c analysis_config.yaml -f $INPUT_FILE -o $OUTPUT_DIR
