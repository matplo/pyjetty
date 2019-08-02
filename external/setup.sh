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
separator "${BASH_SOURCE}"
source ${SCRIPTPATH}/util.sh
setup_python_env

[ "$(get_opt "unset" $@)" == "xyes" ] && unset PYTHONPATH && warning "unsetting PYTHONPATH"
. ${SCRIPTPATH}/setup_lhapdf6.sh 		--version=6.2.3 	 $@
. ${SCRIPTPATH}/setup_hepmc2_cmake.sh 	--version=2.06.09 	 $@
. ${SCRIPTPATH}/setup_hepmc3.sh 		--version=3.1.1  	 $@
. ${SCRIPTPATH}/setup_pythia8.sh 		--version=8235 		 $@
if [ "$(get_opt "install" $@)" == "xyes" ]; then
	note "... running with install"
	. ${SCRIPTPATH}/setup_hepmc3.sh 		--version=3.1.1 --re $@
fi
. ${SCRIPTPATH}/setup_fastjet.sh 		--version=3.3.2 	 $@

for _path in ${HEPMC_DIR} ${HEPMC3_DIR} ${LHAPDF6_DIR} ${PYTHIA8_DIR} ${FASTJET_DIR}
do
	echo ${_path}
	add_path "${_path}/bin"
	_add_python_path=${_path}/lib/python${PYJETTY_PYTHON_VERSION}/site-packages
	if [ "x$(os_darwin)" == "xyes" ]; then
		add_dyldpath "${_path}/lib"
		add_dyldpath "${_path}/lib64"
		add_dyldpath ${_add_python_path}
		add_dyldpath ${_add_python_path}
	else
		add_ldpath "${_path}/lib"
		add_ldpath "${_path}/lib64"
		add_ldpath ${_add_python_path}
		add_ldpath ${_add_python_path}
	fi
	add_pythonpath ${_add_python_path}
done

cd ${cdir}
