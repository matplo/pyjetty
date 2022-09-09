#!/bin/bash
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
	nev=10
fi

./pythia_test_tagger.py --ml  					--py-noue --py-pthatmin 50. --nev ${nev}
./pythia_test_tagger.py --ml --py-hardQCDquarks --py-noue --py-pthatmin 50. --nev ${nev}
./pythia_test_tagger.py --ml --py-hardQCDgluons --py-noue --py-pthatmin 50. --nev ${nev}
