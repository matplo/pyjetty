#! /bin/bash

#SBATCH --job-name="jewel_pp"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=quick
#SBATCH --time=2:00:00
#SBATCH --array=1-1000
#SBATCH --output=/rstorage/ml/egml/data/jewel_pp/slurm-%A_%a.out

JEWEL_DIR=/software/users/james/jewel-2.2.0
JEWEL_EXECUTABLE=jewel-2.2.0-vac
JEWEL_CONFIG_DIR=$JEWEL_DIR/config
JEWEL_CONFIG="params.dat"
OUTDIR=/rstorage/ml/egml/data/jewel_pp/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p $OUTDIR

# Copy config and set output dirs
replace=$OUTDIR/jewel
sed 's!example!'"$replace"'!g' $JEWEL_CONFIG_DIR/$JEWEL_CONFIG >> $OUTDIR/$JEWEL_CONFIG
echo "NJOB $SLURM_ARRAY_TASK_ID" >> $OUTDIR/$JEWEL_CONFIG

echo "Running JEWEL pp"
export LD_LIBRARY_PATH=/software/users/james/lhapdf-5.9.1-install/lib
export LHAPATH=/software/users/james/lhapdf-5.9.1-install/share/lhapdf/PDFsets
$JEWEL_DIR/jewel-2.2.0-vac $OUTDIR/$JEWEL_CONFIG

# Convert hepmc to ntuple
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list
cd /home/james/pyjetty
pipenv run python /software/users/james/pyjetty/pyjetty/alice_analysis/generation/hepmc2antuple_tn_jewel.py -i $OUTDIR/jewel.hepmc -o $OUTDIR/jewel.root --hepmc 2 -g jewel --nev 10000 --no-progress-bar

# Move stdout to appropriate folder
mv /rstorage/ml/egml/data/jewel_pp/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /rstorage/ml/egml/data/jewel_pp/${SLURM_ARRAY_JOB_ID}/
