#! /bin/bash

#SBATCH --job-name=filetest
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-100
#SBATCH --output=/rstorage/alice/AnalysisResults/james/slurm-%A_%a.out

srun process_file_test.sh $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID 100
