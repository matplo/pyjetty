#!/bin/bash

#SBATCH --job-name=djet
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=4:00:00

id
hostname
if [ -d ${1} ]; then
	cd ${1}

	source scl_source enable devtoolset-7
	#module use /home/software/users/ploskon/heppy-devtools7/modules
	#module load heppy/main_python	
	module use /home/software/users/ploskon/pyjetty/modules
	#module load pyjetty/main_python
	module load pyjetty/1.0

	set -x
	/home/software/users/ploskon/heppy/scripts/pipenv_heppy.sh run "${PYJETTY_DIR}/pyjetty/alihfjets/djet.py -f ${2} -o hf_${SLURM_JOB_ID}"
else
	echo "[e] batch job output directory does not exists."
fi
