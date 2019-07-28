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

SCRIPTPATH=$(thisdir)

[ "x${1}" == "xunset" ] && unset PYTHONPATH	&& echo "unsetting PYTHONPATH"
. ${SCRIPTPATH}/setup_hepmc2_cmake.sh
. ${SCRIPTPATH}/setup_lhapdf6.sh
. ${SCRIPTPATH}/setup_pythia8.sh
. ${SCRIPTPATH}/setup_fastjet.sh
. ${SCRIPTPATH}/setup_hepmc3.sh

python_version=$(python3 --version | cut -f 2 -d' ' | cut -f 1-2 -d.)
export PYTHONPATH=${FASTJET_DIR}/lib/python${python_version}/site-packages:${PYTHONPATH}
export PYTHONPATH=${HEPMC2_DIR}/lib:${PYTHONPATH}
export PYTHONPATH=${LHAPDF6_DIR}/lib/python${python_version}/site-packages:${PYTHONPATH}
export PYTHONPATH=${PYTHIA_DIR}/lib:${PYTHONPATH}

export PATH=${HEPMC_DIR}/bin:${LHAPDF6_DIR}/bin:${PYTHIA8_DIR}/bin:${FASTJET_DIR}/bin:${PATH}
if [ -z ${LD_LIBRARY_PATH} ]; then
	export LD_LIBRARY_PATH=${HEPMC_DIR}/lib:${LHAPDF6_DIR}/lib:${PYTHIA_DIR}/lib:${FASTJET_DIR}/lib
else
	export LD_LIBRARY_PATH=${HEPMC_DIR}/lib:${LHAPDF6_DIR}/lib:${PYTHIA_DIR}/lib:${FASTJET_DIR}/lib:${LD_LIBRARY_PATH}
fi
if [ -z ${DYLD_LIBRARY_PATH} ]; then
	export DYLD_LIBRARY_PATH=${HEPMC_DIR}/lib:${LHAPDF6_DIR}/lib:${PYTHIA_DIR}/lib:${FASTJET_DIR}/lib
else
	export DYLD_LIBRARY_PATH=${HEPMC_DIR}/lib:${LHAPDF6_DIR}/lib:${PYTHIA_DIR}/lib:${FASTJET_DIR}/lib:${DYLD_LIBRARY_PATH}
fi
