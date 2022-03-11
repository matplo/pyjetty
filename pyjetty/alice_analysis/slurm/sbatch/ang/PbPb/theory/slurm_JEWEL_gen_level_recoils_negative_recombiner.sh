#! /bin/bash

#SBATCH --job-name=jwl_rec_nosub
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --partition=quick
#SBATCH --time=2:00:00
#SBATCH --array=1-400
#SBATCH --output=/rstorage/alice/AnalysisResults/ang/slurm-%A_%a.out

FILE_PATHS='/rstorage/generators/jewel_alice/tree_gen/851894/files.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 7 nodes * 20 cores active
FILES_PER_JOB=20
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
  srun process_JEWEL_gen_level_recoils_negative_recombiner.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done

