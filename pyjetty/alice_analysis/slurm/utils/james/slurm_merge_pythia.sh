#! /bin/bash

#SBATCH --job-name=mergepthat
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=10
#SBATCH --partition=quick
#SBATCH --time=2:00:00
#SBATCH --array=1-20
#SBATCH --output=/rstorage/alice/AnalysisResults/james/slurm-%A_%a.out

srun merge_pythia.sh $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
