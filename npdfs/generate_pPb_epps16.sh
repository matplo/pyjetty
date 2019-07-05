#!/bin/bash

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
	outdir="${PWD}/epps16-pPb-ehigh-bias-${biaspower}-${biasref}"
	mkdir -p ${outdir}
	sets=$(seq -1 0)
	parallel --bar ./pPb_epps16.py -g ${outdir}/pPb_npdf_compare.dat --epps16set {} --ecm high --biaspow 4 --biasref 20 -n 1000 ::: ${sets}
	parallel --bar ./pPb_epps16.py -g ${outdir}/pPb_npdf_compare_charm.dat --charm --epps16set {} --ecm high --biaspow 4 --biasref 20 -n 1000 ::: ${sets}
	parallel --bar ./pPb_epps16.py -g ${outdir}/pPb_npdf_compare_photon.dat --photon --epps16set {} --ecm high --biaspow 4 --biasref 20 -n 1000 ::: ${sets}
else
	pthatmin=${1}
	if [ ! -z ${pthatmin} ]; then
		outdir=${PWD}/epps16-pPb-ehigh-pthat${pthatmin}/
		mkdir -p ${outdir}
		sets=$(seq -1 0)
		parallel --bar ./pPb_epps16.py -g ${outdir}/pPb_npdf_compare.dat --epps16set {} --ecm high --pthatmin ${pthatmin} -n 1000 ::: ${sets}
		parallel --bar ./pPb_epps16.py -g ${outdir}/pPb_npdf_compare_charm.dat --charm --epps16set {} --ecm high --pthatmin ${pthatmin} -n 1000 ::: ${sets}
		parallel --bar ./pPb_epps16.py -g ${outdir}/pPb_npdf_compare_photon.dat --photon --epps16set {} --ecm high --pthatmin ${pthatmin} -n 1000 ::: ${sets}
	fi
fi
