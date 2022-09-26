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

for seed in 123456
do
	# quarks and gluons
	${THISD}/pythia_parton_tag.py --ml 						          --py-noue --py-pthatmin 500 --jet-ptmin 500 --jet-R 0.8 --jet-ptmax 550 $@
	# only ff->gg
	${THISD}/pythia_parton_tag.py --ml --py-hardQCDgluons 	--py-noue --py-pthatmin 500 --jet-ptmin 500 --jet-R 0.8 --jet-ptmax 550 $@
	# only ff->qq or ff->qqbar
	${THISD}/pythia_parton_tag.py --ml --py-hardQCDquarks 	--py-noue --py-pthatmin 500 --jet-ptmin 500 --jet-R 0.8 --jet-ptmax 550 $@
done