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
[ "x$(get_opt "python2" $@)" == "xyes" ] && export PYJETTY_USER_PYTHON_VERSION=python2
[ "x$(get_opt "python3" $@)" == "xyes" ] && export PYJETTY_USER_PYTHON_VERSION=python3
setup_python_env

[ "x$(get_opt "unset" $@)" == "xyes" ] && unset PYTHONPATH && warning "unsetting PYTHONPATH"
. ${SCRIPTPATH}/setup_lhapdf6.sh 		--version=6.2.3 	 $@
. ${SCRIPTPATH}/setup_hepmc2_cmake.sh 	--version=2.06.09 	 $@
. ${SCRIPTPATH}/setup_hepmc3.sh 		--version=3.1.1  	 $@
. ${SCRIPTPATH}/setup_pythia8.sh 		--version=8235 		 $@
if [ "$(get_opt "install" $@)" == "xyes" ]; then
	note "... running with install"
	. ${SCRIPTPATH}/setup_hepmc3.sh 		--version=3.1.1 --re $@
fi
. ${SCRIPTPATH}/setup_fastjet.sh 		--version=3.3.2 	 $@

for path in ${HEPMC2_DIR} ${HEPMC3_DIR} ${LHAPDF6_DIR} ${PYTHIA8_DIR} ${FASTJET_DIR}
do
	note "using ${path} to setup environment..."
	bin_path="${path}/bin"
	lib_path="${path}/lib"
	lib64_path="${path}/lib64"
	python_path="${path}/lib/python${PYJETTY_PYTHON_VERSION}/site-packages"
	python_path64="${path}/lib64/python${PYJETTY_PYTHON_VERSION}/site-packages"
	add_path ${bin_path}
	if [ "x$(os_darwin)" == "xyes" ]; then
		add_dyldpath "${lib_path}"
		add_dyldpath "${lib64_path}"
		add_dyldpath "${python_path}"
		add_dyldpath "${python_path64}"
	else
		add_ldpath "${lib_path}"
		add_ldpath "${lib64_path}"
		add_ldpath "${python_path}"
		add_ldpath "${python_path64}"
	fi
	add_pythonpath "${python_path}"
	add_pythonpath "${python_path64}"
	add_pythonpath "${lib_path}"
	add_pythonpath "${lib64_path}"
done

cd ${cdir}
