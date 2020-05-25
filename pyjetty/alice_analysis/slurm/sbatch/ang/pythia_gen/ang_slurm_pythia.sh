#! /bin/bash

#SBATCH --job-name="pythiagen"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-280
#SBATCH --output=/rstorage/alice/AnalysisResults/ang/slurm-%A_%a.out

# Center of mass energy in GeV
ECM=5020

# Number of events per pT-hat bin (for statistics)
NEV_DESIRED=1400000

# Array of the jet radii for which to generate events
JETR_ARR=(0.2 0.4)

# Lower edges of the pT-hat bins
PTHAT_BINS=(5 7 9 12 16 21 28 36 45 57 70 85 99 115 132 150 169 190 212 235)
echo "Number of pT-hat bins: ${#PTHAT_BINS[@]}"

# Currently we have 8 nodes * 20 cores active
NCORES=280
NEV_PER_JOB=$(( $NEV_DESIRED * ${#PTHAT_BINS[@]} / $NCORES ))
echo "Number of events per job: $NEV_PER_JOB"
NCORES_PER_BIN=$(( $NCORES / ${#PTHAT_BINS[@]} ))
echo "Number of cores per pT-hat bin: $NCORES_PER_BIN"

BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) / $NCORES_PER_BIN + 1))
CORE_IN_BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) % $NCORES_PER_BIN + 1))
PTHAT_MIN=${PTHAT_BINS[$(( $BIN - 1 ))]}
if [ $BIN -lt ${#PTHAT_BINS[@]} ]; then
    USE_PTHAT_MAX=true
	PTHAT_MAX=${PTHAT_BINS[$BIN]}
	echo "Calculating bin $BIN (pThat=[$PTHAT_MIN,$PTHAT_MAX]) with core number $CORE_IN_BIN"
else
    USE_PTHAT_MAX=false
	echo "Calculating bin $BIN (pThat_min=$PTHAT_MIN) with core number $CORE_IN_BIN"
fi

SEED=$(( ($CORE_IN_BIN - 1) * NEV_PER_JOB + 1111 ))

# Do the PYTHIA simulation & matching
OUTDIR="/rstorage/alice/AnalysisResults/ang/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN"
mkdir -p $OUTDIR
cd /home/ezra/pyjetty/pyjetty/alice_analysis/

for JETR in "${JETR_ARR[@]}"
do
	# First do with MPI off
	FNAME="PythiaResults_R${JETR}_MPIoff.root"
	if $USE_PTHAT_MAX; then
		python process/user/ang_pp/pythia_parton_hadron.py --output "$OUTDIR/$FNAME" \
			--user-seed $SEED --jetR $JETR --py-pthatmin $PTHAT_MIN --py-ecm $ECM --py-noMPI \
			--nev $NEV_PER_JOB --pythiaopts PhaseSpace:pTHatMax=$PTHAT_MAX
	else
		python process/user/ang_pp/pythia_parton_hadron.py --output "$OUTDIR/$FNAME" \
			--user-seed $SEED --jetR $JETR --py-pthatmin $PTHAT_MIN --py-ecm $ECM \
			--py-noMPI --nev $NEV_PER_JOB
	fi

	# ... then with MPI on
	FNAME="PythiaResults_R${JETR}.root"
	if $USE_PTHAT_MAX; then
		python process/user/ang_pp/pythia_parton_hadron.py --output "$OUTDIR/$FNAME" \
			--user-seed $SEED --jetR $JETR --py-pthatmin $PTHAT_MIN --py-ecm $ECM \
			--nev $NEV_PER_JOB --pythiaopts PhaseSpace:pTHatMax=$PTHAT_MAX
	else
		python process/user/ang_pp/pythia_parton_hadron.py --output "$OUTDIR/$FNAME" \
			--user-seed $SEED --jetR $JETR --py-pthatmin $PTHAT_MIN \
			--py-ecm $ECM --nev $NEV_PER_JOB
	fi
done
