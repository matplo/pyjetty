#! /bin/bash

#SBATCH --job-name="Rg_herwig"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=3
#SBATCH --partition=long
#SBATCH --time=48:00:00
#SBATCH --array=1-2560
#SBATCH --output=/rstorage/alice/AnalysisResults/theta_g/slurm-%A_%a.out

FILE_PATHS='/rstorage/alice/sim/herwig_gen/files.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=$(( $NFILES / 2560 + 1 ))
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
  srun theta_g_gen.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
