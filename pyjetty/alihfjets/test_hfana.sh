#!/bin/bash -l

#SBATCH --job-name=test_hfana
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00

id
hostname
date

source scl_source enable devtoolset-7
module use /home/software/users/${USER}/pyjetty/modules
module avail
module load pyjetty/1.0
module list

set -x
executable=${1}
flist=${2}
outputdir=${3}
output=${outputdir}/$(basename ${flist})

if [ -e ${flist} ]; then
	mkdir -p ${outputdir}
	if [ -d ${outputdir} ]; then
		cd ${outputdir}
		if [ -x ${executable} ]; then
			${HEPPY_DIR}/scripts/pipenv_heppy.sh run ${executable} -f ${flist} -o ${output}
		else
			echo "[e] no executable: ${executable}"
		fi
	else
		echo "[e] no outputdir: ${outputdir}"
	fi
else
	echo "[e] no file list: ${flist}"
fi

date