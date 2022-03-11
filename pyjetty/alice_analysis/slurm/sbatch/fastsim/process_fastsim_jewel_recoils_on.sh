#! /bin/bash

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
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f6-8)
echo $INPUT_FILE
echo $OUTPUT_SUFFIX
OUTPUT_DIR="/rstorage/generators/jewel_alice/tree_fastsim_recoils_on/$JOB_ID/$OUTPUT_SUFFIX"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
module use /home/ezra/heppy/modules
module load heppy/1.0
module use /home/ezra/pyjetty/modules
module load pyjetty/1.0
module list

# Run main script
cd /home/ezra/pyjetty/pyjetty/alice_analysis/process/user/fastsim
python eff_smear.py -i $INPUT_FILE -o $OUTPUT_DIR --jewel

# Move stdout to appropriate folder
mv /rstorage/generators/jewel_alice/tree_fastsim_recoils_on/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/generators/jewel_alice/tree_fastsim_recoils_on/${JOB_ID}/
