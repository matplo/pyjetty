#!/usr/bin/env bash

function thisdir()
{
        SOURCE="${BASH_SOURCE[0]}"
        while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
          DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
          SOURCE="$(readlink "$SOURCE")"
          [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
        done
        DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
        echo ${DIR}
}
THISD=$(thisdir)

source ${HEPPY_DIR}/scripts/util.sh

nev=$(get_opt "nev" $@)
if [ -z ${nev} ]; then
	nev=1000
fi

for seed in 123456
do
  # quarks and gluons
  ${THISD}/pythia_gen.py                     --py-noue --py-pthatmin 100. --py-seed ${seed} --ml --hepmc --nev ${nev}
  # only ff->gg
	${THISD}/pythia_gen.py --py-hardQCDgluons --py-noue --py-pthatmin 100. --py-seed ${seed} --ml --hepmc --nev ${nev}
  # only ff->qq or ff->qqbar
	${THISD}/pythia_gen.py --py-hardQCDquarks --py-noue --py-pthatmin 100. --py-seed ${seed} --ml --hepmc --nev ${nev}
done
