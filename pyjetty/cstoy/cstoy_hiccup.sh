#!/bin/bash

#SBATCH --job-name=cstoy
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=4:00:00

id
hostname
mkdir -p /storage/u/ploskon/cstoy
cd /storage/u/ploskon/cstoy

source scl_source enable devtoolset-7
#module use /home/software/users/ploskon/heppy-devtools7/modules
#module load heppy/main_python	
module use /home/software/users/ploskon/pyjetty/modules
#module load pyjetty/main_python
module load pyjetty/1.0

set -x
/home/software/users/ploskon/heppy/scripts/pipenv_heppy.sh run "$PYJETTYDIR/pyjetty/cstoy/cstoy.py --nev ${3} --noue --pthatmin 100 --alpha ${1} --dRmax ${2}"
