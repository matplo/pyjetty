#! /bin/bash

#SBATCH --job-name=ml-compute
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=quick
#SBATCH --time=2:00:00
#SBATCH --array=1-100
#SBATCH --output=/rstorage/ml/egml/nsubjettiness/slurm-%A_%a.out

FILE_PATHS='/rstorage/ml/egml/data/pyquen/1_0_1_100k_0fm/files.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=1
echo "Files per job: $FILES_PER_JOB"

STOP=$(( SLURM_ARRAY_TASK_ID*FILES_PER_JOB ))
START=$(( $STOP - $(( $FILES_PER_JOB - 1 )) ))

if (( $STOP > $NFILES ))
then
  STOP=$NFILES
fi

echo "START=$START"
echo "STOP=$STOP"

for (( JOB_N = $START; JOB_N <= $STOP; JOB_N++ ))
do
  FILE=$(sed -n "$JOB_N"p $FILE_PATHS)
  srun process_compute_nsubjettiness_pyquen.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
