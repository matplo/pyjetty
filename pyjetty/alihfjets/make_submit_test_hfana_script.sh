#!/bin/bash

cpwd=${PWD}

pp5TeVfiles=/rstorage/alice/data/LHC17pq/filelist.txt
flist=${pp5TeVfiles}

dname=$(date +"%Y-%m-%d-%H-%M")

flistdname=$(dirname ${flist})
outputdir=$(basename ${flistdname})
outputdir=/rstorage/${USER}/test_hfana/${outputdir}/${dname}
mkdir -pv ${outputdir}
cd ${outputdir}
pwd

cp -v ${flist} .
split --additional-suffix=.flist -d -l 5 -a 5 ${flist}

job_lists=$(find $PWD -name "*.flist")

cp -v ${PYJETTY_DIR}/pyjetty/alihfjets/test_hfana.sh .
cp -v ${PYJETTY_DIR}/pyjetty/alihfjets/test_hfana.py .

executable=${PWD}/test_hfana.py

submit_script=${PWD}/submit_all.sh
rm -f ${submit_script}

for fl in ${job_lists}
do
	echo "sbatch --chdir=${PWD} --output=${fl}.output --error=${fl}.error ${PWD}/test_hfana.sh ${executable} ${fl} ${outputdir}" | tee -a ${submit_script}
	chmod +x ${submit_script}
done

echo "[i] created: ${submit_script}"

cd ${cpwd}
