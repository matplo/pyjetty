#!/bin/bash

fname=${1}

outdir_date=$(date +"%Y-%m-%d-at-%H.%M.%S")
outdir=/rstorage/ploskon/alihfjets/${outdir_date}

if [ ! -z ${1} ] && [ -e ${1} ]; then
	mkdir -p ${outdir}
	file_list=$(cat ${1})
	for fname in ${file_list}
	do
		if [ -e ${fname} ]; then
			echo ${fname}
			sbatch -D ${outdir} ./djet_hiccup.sh ${outdir} ${fname}
		fi
	done
fi