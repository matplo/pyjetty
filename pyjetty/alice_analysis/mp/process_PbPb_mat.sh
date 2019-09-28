#! /bin/bash

#SBATCH --job-name=matPbPb
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=4:00:00

# Take two command line arguments -- (1) input file path, (2) output dir prefix
if [ ! -z "$1" ]; then
  INPUT_FILE=$1
  #echo "Input file: $INPUT_FILE"
else
  echo "Wrong command line arguments"
fi

if [ ! -z "$2" ]; then
  OUTPUT_PREFIX=$2
else 
  echo "Wrong command line arguments"
fi

# Load modules
source scl_source enable devtoolset-7
module use /home/software/users/ploskon/heppy-devtools7/modules
module load heppy/main_python
module use /home/software/users/ploskon/pyjetty/modules
module load pyjetty/main_python

# source /home/software/users/ploskon/RooUnfold-gitlab/build/setup.sh
# export ROOUNFOLDDIR=/home/software/users/ploskon/RooUnfold-gitlab/build

# Run python script via pipenv
##cd /software/users/james/pyjetty/pyjetty/alice_rg
##pipenv run python process_rg_mc.py -c analysis_config.yaml -f $INPUT_FILE -b $BIN -w $SCALE_FACTOR_PATH -o $OUTPUT_DIR
${HEPPYDIR}/scripts/pipenv_heppy.sh run "${PYJETTYDIR}/pyjetty/alice_analysis/mp/matPbPb.py --input ${INPUT_FILE} --outprefix ${OUTPUT_PREFIX}"
