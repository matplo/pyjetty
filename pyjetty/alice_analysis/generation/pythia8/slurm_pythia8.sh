#! /bin/bash

#SBATCH --job-name="pythia8"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=quick
#SBATCH --time=2:00:00
#SBATCH --array=1-1000
#SBATCH --output=/rstorage/ml/egml/data/pythia8/slurm-%A_%a.out

BASE_DIR=/software/users/james/pyjetty/pyjetty/alice_analysis/generation
PYTHIA_DIR=$BASE_DIR/pythia8
PYTHIA_CONFIG=$PYTHIA_DIR/settings.cmnd
OUTDIR=/rstorage/ml/egml/data/pythia8/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p $OUTDIR

# Enter environment
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list

# Run PYTHIA
echo "Running PYTHIA8"
cd /home/james/pyjetty
pipenv run python $PYTHIA_DIR/pythia_gen_write_hepmc.py --py-cmnd $PYTHIA_CONFIG --nev 1000 -o $OUTDIR

# Convert hepmc to ntuple
pipenv run python $BASE_DIR/hepmc2antuple_tn.py -i $OUTDIR/pythia8.hepmc -o $OUTDIR/pythia8.root --hepmc 2 -g pythia --nev 1000 --no-progress-bar

# Move stdout to appropriate folder
mv /rstorage/ml/egml/data/pythia8/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /rstorage/ml/egml/data/pythia8/${SLURM_ARRAY_JOB_ID}/
