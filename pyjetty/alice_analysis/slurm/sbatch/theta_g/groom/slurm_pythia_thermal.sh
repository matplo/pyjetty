#! /bin/bash

#SBATCH --job-name=pythia-emb
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-1000
#SBATCH --output=/rstorage/james/groom/AnalysisResults/slurm-%A_%a.out

# Note: Set a separate job for each file -- otherwise they will rewrite each other

FILE_PATHS='/rstorage/ploskon/groom/pythia8/pppthatbias-root/files_30XXX.txt'
#FILE_PATHS='/rstorage/ploskon/groom/pythia8/pppthatbias-root/files_31XXX.txt'
#FILE_PATHS='/rstorage/ploskon/groom/pythia8/pppthatbias-root/files_32XXX.txt'
#FILE_PATHS='/rstorage/ploskon/groom/pythia8/pthatminX-root/files_50.txt'
#FILE_PATHS='/rstorage/ploskon/groom/pythia8/pthatminX-root/files_100.txt'
#FILE_PATHS='/rstorage/ploskon/groom/pythia8/pthatminX-root/files_200.txt'

NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 7 nodes * 20 cores active
#FILES_PER_JOB=$(( $NFILES / 10 ))
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
  srun process_pythia_thermal.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
