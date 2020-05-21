#! /bin/bash

#SBATCH --job-name=james-convert
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-10
#SBATCH --output=/rstorage/james/groom/slurm-%A_%a.out

FILE_PATHS='/rstorage/ploskon/groom/pythia8/pthatmin50/hepmc_filelist.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 7 nodes * 20 cores active
FILES_PER_JOB=$(( $NFILES / 10 ))
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
  srun process_convert_pythia_hepmc.sh $FILE
done
