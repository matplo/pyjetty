#! /bin/bash

#SBATCH --job-name="JEWEL-ALICE"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --time=24:00:00
#SBATCH --array=1-8000
#SBATCH --output=/rstorage/generators/jewel_alice/slurm-%A_%a.out

# Number of events per pT-hat bin
NEV_DESIRED=2000000

# Lower edges of the pT-hat bins
PTHAT_BINS=(5 7 9 12 16 21 28 36 45 57 70 85 99 115 132 150 169 190 212 235)
echo "Number of pT-hat bins: ${#PTHAT_BINS[@]}"

# Set number of jobs to split generation between, and calculate how many jobs per pt-hat bin
NCORES=8000
NEV_PER_JOB=$(( $NEV_DESIRED * ${#PTHAT_BINS[@]} / $NCORES ))
echo "Number of events per job: $NEV_PER_JOB"
NCORES_PER_BIN=$(( $NCORES / ${#PTHAT_BINS[@]} ))
echo "Number of cores per pT-hat bin: $NCORES_PER_BIN"

BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) / $NCORES_PER_BIN + 1))
CORE_IN_BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) % $NCORES_PER_BIN + 1))

# Set JEWEL paths, and create output directory
JEWEL_DIR=/software/users/james/jewel-2.2.0
JEWEL_CONFIG_DIR=$JEWEL_DIR/config
JEWEL_EXECUTABLE=jewel-2.2.0-simple
OUTDIR_BASE="/rstorage/generators/jewel_alice"
OUTDIR=$OUTDIR_BASE/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN
mkdir -p $OUTDIR

# Construct config file for this pt-hat bin -- set log,hepmc output locations, and set n_event and pt-hat bin
replace=$OUTDIR/jewel
sed 's!example!'"$replace"'!g' $JEWEL_CONFIG_DIR/params.dat >> $OUTDIR/params.dat

PTHAT_MIN=${PTHAT_BINS[$(( $BIN - 1 ))]}
echo $'\n' >> $OUTDIR/params.dat
echo "NEVENT $NEV_PER_JOB" >> $OUTDIR/params.dat
echo "PTMIN $PTHAT_MIN" >> $OUTDIR/params.dat
if [ $BIN -lt ${#PTHAT_BINS[@]} ]; then
	PTHAT_MAX=${PTHAT_BINS[$BIN]}
    echo "PTMAX $PTHAT_MAX" >> $OUTDIR/params.dat
	echo "Calculating bin $BIN (pThat=[$PTHAT_MIN,$PTHAT_MAX]) with core number $CORE_IN_BIN"
else
	echo "Calculating bin $BIN (pThat_min=$PTHAT_MIN) with core number $CORE_IN_BIN"
fi
echo "NJOB $SLURM_ARRAY_TASK_ID" >> $OUTDIR/params.dat

echo "Running JEWEL PbPb"
export LD_LIBRARY_PATH=/software/users/james/lhapdf-5.9.1-install/lib
export LHAPATH=/software/users/james/lhapdf-5.9.1-install/share/lhapdf/PDFsets
$JEWEL_DIR/$JEWEL_EXECUTABLE $OUTDIR/params.dat

# Convert hepmc to ntuple
module use /software/users/james/heppy/modules
module load heppy/1.0
module use /software/users/james/pyjetty/modules
module load pyjetty/1.0
module list
cd /home/james/pyjetty
pipenv run python /software/users/james/pyjetty/pyjetty/alice_analysis/generation/hepmc2antuple_tn_jewel.py -i $OUTDIR/jewel.hepmc -o $OUTDIR/jewel.root --hepmc 2 -g jewel_charged --nev $NEV_PER_JOB --no-progress-bar

# Move stdout to appropriate folder
mv $OUTDIR_BASE/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $OUTDIR_BASE/${SLURM_ARRAY_JOB_ID}/