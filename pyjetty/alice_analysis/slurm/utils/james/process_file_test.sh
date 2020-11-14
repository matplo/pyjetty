#! /bin/bash

#OUTPUTDIR=/rstorage/alice/data/LHC18qr/413-414
#OUTPUTDIR=/rstorage/alice/data/LHC19f4/403
#OUTPUTDIR=/rstorage/alice/data/LHC17pq/448
OUTPUTDIR=/rstorage/alice/data/LHC18b8/520
#OUTPUTDIR=/rstorage/alice/data/LHC18bdeghijkklnop/436

FILELIST=$OUTPUTDIR/files.txt

if [ "$1" != "" ]; then
  JOB_ID=$1
  echo "Job ID: $JOB_ID"
else 
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  TASK_ID=$2
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

if [ "$3" != "" ]; then
  NJOBS=$3
  echo "NJOBS: $NJOBS"
else
  echo "Wrong command line arguments"
fi

# Load modules
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list

# Run python script via pipenv
cd /home/james/pyjetty
pipenv run python pyjetty/alice_analysis/slurm/utils/james/test_ROOT_trees.py -f $FILELIST -n $NJOBS -i $TASK_ID -o $OUTPUTDIR
