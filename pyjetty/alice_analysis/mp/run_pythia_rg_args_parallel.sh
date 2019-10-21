#!/bin/bash

isr="--py-noISR"
mpi="--py-noMPI"
had="--pythiaopts HadronLevel:all=off"
pthatmin="20 40 60 80"
jetR="0.4"

pthats="4 11 21 36 56 84 117 156 200 249 1000"
pthats="5 7 9 12 16 21 28 36 45 57 70 85 99 115 132 150 169 190 212 235"

run=${1}
overw=""
[ "x${2}" == "xoverwrite" ] && overw="--overwrite"
if [ "x${run}" == "xrun" ]; then
	njobs=${2}
	dryrun=""
	[ "x${2}" == "xdry" ] && dryrun="--dry-run" && njobs=${3}
	[ "x${3}" == "xdry" ] && dryrun="--dry-run"
	[ -z ${njobs} ] && njobs=1
	runs=$(seq 1 ${njobs})
	nev=25000
	cmnd="${PYJETTYDIR}/pyjetty/alice_analysis/mp/pythia_rg.py --nev ${nev} --jetR ${jetR} --py-ecm 5020 --charged --RreclusterR0"
	parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1} 	                 --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1} ${had}               --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1} ${had} ${mpi}        --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1} ${had}        ${isr} --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1} ${had} ${mpi} ${isr} --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1}        ${mpi}        --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1}        ${mpi} ${isr} --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
	#parallel ${dryrun} --joblog hjet_simple_hTT_softQCD.log --keep-order --tag "${cmnd} --runid {1}               ${isr} --py-pthatmin {2}" ::: ${runs} ::: ${pthats}
fi


