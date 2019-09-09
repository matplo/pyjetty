#! /bin/bash

#SBATCH --job-name=rgtest
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=quick
#SBATCH --time=1:00:00
#SBATCH --array=1-5
#SBATCH --output=/storage/u/alice/AnalysisResults/slurm-%A_%a.out

#5291

FILE_PATH='/rstorage/u/alice/LHC17pq/145/files_test.txt'
NFILES=$(wc -l $FILE_PATH)
echo "N files: ${NFILES}"

FILES=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_PATH)

srun process_rg_data.sh $FILES
