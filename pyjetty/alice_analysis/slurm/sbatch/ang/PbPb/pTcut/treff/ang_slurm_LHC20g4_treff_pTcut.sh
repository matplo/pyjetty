#! /bin/bash

#SBATCH --job-name="ang_treff"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=3
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-800
#SBATCH --output=/rstorage/alice/AnalysisResults/ang/slurm-%A_%a.out

FILE_PATHS='/rstorage/alice/data/LHC20g4/568/files.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=$(( $NFILES / 800 + 1 ))
echo "Files per job: $FILES_PER_JOB"

STOP=$(( SLURM_ARRAY_TASK_ID * FILES_PER_JOB ))
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
  srun ang_LHC20g4_treff_pTcut.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
