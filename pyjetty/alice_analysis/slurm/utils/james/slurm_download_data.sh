#! /bin/bash

#SBATCH --job-name=download
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=10
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-14
#SBATCH --output=/rstorage/alice/AnalysisResults/james/slurm-%A_%a.out

srun downloadData.sh $SLURM_ARRAY_TASK_ID
sleep 1
