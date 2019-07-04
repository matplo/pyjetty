#!/bin/bash

function run_reco_jets()
{
	files=${1}
	nev=
	for fn in ${files}
	do
		echo "[${fn}]"
		if [ -z ${nev} ]; then
			nev=$(cat ${fn} | grep "E " | wc -l)
		fi
		# ./reco_jets_epps16.py --read ${fn} -n ${nev}	
	done
}
export -f run_reco_jets

indir=${1}
nev=${2}
if [ ! -z ${indir} ]; then
	sfiles=$(find ${indir} -name "*.dat")
	# find ${indir} -name "*.dat" | parallel -X --dryrun run_something
	# find ${indir} -name "*.dat" | parallel -X ./reco_jets_epps16.py --read {}
	# echo ${sfiles} | parallel -X run_reco_jets {}

	# parallel --bar ./reco_jets_epps16.py --read {} -n $(cat ${fn} | grep "E " | wc -l) ::: ${sfiles}
	fn=$(echo ${sfiles}  | head -n1 | cut -d " " -f1)
	if [ -z ${nev} ]; then
		nev=$(cat ${fn} | grep "E " | wc -l)
	fi
	echo "number of events per file: ${nev}"
	parallel --bar ./reco_jets_epps16.py --read {} -n ${nev}  ::: ${sfiles}
else
	echo "usage: ${0} <input_directory> [nev]"
fi

#indir=${1}
#nev=${2}
## sfiles=$(find ./epps16-pp-5000-bias -name "*.dat")
## sfiles=$(find ./epps16-pp-5000-pthat5 -name "*.dat")
#if [ ! -z ${indir} ]; then
#	sfiles=$(find ${indir} -name "*.dat")
#	for fn in ${sfiles}
#	do
#		echo "[i] processing ${fn}"
#		if [ -z ${nev} ]; then
#			nev=$(cat ${fn} | grep "E " | wc -l)
#		fi
#		./reco_jets_epps16.py --read ${fn} -n ${nev}
#	done
#else
#	echo "usage: $(basename ${0}) <directory> [nevents]"
#fi
