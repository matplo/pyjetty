#!/bin/bash

cdir=$(pwd)

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
source ${SCRIPTPATH}/util.sh
separator "${BASH_SOURCE}"

[ "$(get_opt "unset" $@)" == "xyes" ] && unset PYTHONPATH && warning "unsetting PYTHONPATH"
. ${SCRIPTPATH}/setup_lhapdf6.sh 		--version=6.2.3 	 $@
. ${SCRIPTPATH}/setup_hepmc2_cmake.sh 	--version=2.06.09 	 $@
. ${SCRIPTPATH}/setup_hepmc3.sh 		--version=3.1.1  	 $@
. ${SCRIPTPATH}/setup_pythia8.sh 		--version=8235 		 $@
if [ "$(get_opt "install" $@)" == "xyes" ]; then
	warning "... running with install"
	. ${SCRIPTPATH}/setup_hepmc3.sh 		--version=3.1.1 --re $@
fi
. ${SCRIPTPATH}/setup_fastjet.sh 		--version=3.3.2 	 $@

python_version=$(python3 --version | cut -f 2 -d' ' | cut -f 1-2 -d.)
export PYTHONPATH=${FASTJET_DIR}/lib/python${python_version}/site-packages:${PYTHONPATH}
# export PYTHONPATH=${HEPMC2_DIR}/lib:${PYTHONPATH}
export PYTHONPATH=${HEPMC3_DIR}/lib:${PYTHONPATH}
export PYTHONPATH=${LHAPDF6_DIR}/lib/python${python_version}/site-packages:${PYTHONPATH}
export PYTHONPATH=${PYTHIA_DIR}/lib:${PYTHONPATH}

for _path in ${HEPMC_DIR} ${HEPMC3_DIR} ${LHAPDF6_DIR} ${PYTHIA8_DIR} ${FASTJET_DIR}
do
	echo ${_path}
	if [ ! -z ${_path} ] && [ -d ${_path} ]; then
		echo_info "adding ${_path}"
		if [ -z ${PATH} ]; then
			export PATH=${_path}/bin
		else
			export PATH=${_path}/bin:${PATH}
		fi
		if [ -z ${LD_LIBRARY_PATH} ]; then
			export LD_LIBRARY_PATH=${_path}/lib
		else
			export LD_LIBRARY_PATH=${_path}/lib:${LD_LIBRARY_PATH}
		fi
		if [ -z ${DYLD_LIBRARY_PATH} ]; then
			export DYLD_LIBRARY_PATH=${_path}/lib
		else
			export DYLD_LIBRARY_PATH=${_path}/lib:${DYLD_LIBRARY_PATH}
		fi
	fi
done

cd ${cdir}