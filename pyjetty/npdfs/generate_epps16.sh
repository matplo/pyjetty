#!/bin/bash

# outdir=${PWD}/epps16-pp-5000-bias/
# mkdir -p ${outdir}
# ./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --allset --ecm 5000. --biaspow 4 --biasref 20

# outdir=${PWD}/epps16-pp-5000-pthat5/
# mkdir -p ${outdir}
# ./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --allset --ecm 5000. --pthatmin 5 -n 10000

smode=${1}

if [ -z ${smode} ]; then
	echo $(basename ${0}) "[<pthatmin> | bias] [biaspower] [biasref]"
	exit 0
fi

if [ "xbias" == "x$smode" ]; then
	biaspower=4
	[ ! -z ${2} ] && biaspower=${2}
	biasref=20
	[ ! -z ${3} ] && biasref=${3}
	echo "mode : ${smode} && biaspower : ${biaspower} && biasref : ${biasref}"
	outdir="${PWD}/epps16-pp-5000-bias-${biaspower}-${biasref}"
	mkdir -p ${outdir}
	# for i in $(seq -1 40)
	# do
	# 	echo "set ${i}"
	# 	./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --epps16set ${i} --ecm 5000. --biaspow 4 --biasref 20 -n 10000
	# done
	sets=$(seq -1 40)
	parallel --bar ./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --epps16set {} --ecm 5000. --biaspow 4 --biasref 20 -n 10000 ::: ${sets}
else
	pthatmin=${2}
	if [ ! -z ${pthatmin} ]; then
		outdir=${PWD}/epps16-pp-5000-pthat${pthatmin}/
		mkdir -p ${outdir}
		for i in $(seq -1 40)
		do
			echo "set ${i}"
			./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --epps16set ${i} --ecm 5000. --pthatmin ${pthatmin} -n 10000
		done
		# ./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --allset --ecm 5000. --pthatmin ${pthatmin} -n 10000
	fi
fi