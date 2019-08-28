#!/bin/bash

source ${HEPPYDIR}/scripts/util.sh
nev=$(get_opt "nev" $@)

outputfiles=""
for pthat in 50 100 500
do
	outputfile=pythia_hardQCD_${pthat}.dat
	[ ! -f ${outputfile} ] && ./pythia_gen_write_hepmc.py --pthatmin ${pthat} --hardQCD --output ${outputfile} --nev ${nev}
done

outputfile=pythia_hardQCD_w.dat
[ ! -f ${outputfile} ] && ./pythia_gen_write_hepmc.py --biaspow 4 --biasref 50 --hardQCD --output ${outputfile} --nev ${nev}

reco=$(get_opt "reco" $@)
if [ ${reco} ]; then
	./hepmc_jewel_jetreco_lundjets_joblib.py -d . -p "*.dat" --nev ${nev} --clean
fi
