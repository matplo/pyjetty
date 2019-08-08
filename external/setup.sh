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
THISD=$(thisdir)
source ${THISD}/../scripts/util.sh
separator "${BASH_SOURCE}"
build_with_python="python"
[ "x$(get_opt "python2" $@)" == "xyes" ] && build_with_python="python2"
[ "x$(get_opt "python3" $@)" == "xyes" ] && build_with_python="python3"
${THISD}/../scripts/make_module.sh ${build_with_python}
module use ${THISD}/../modules
module avail
separator ''
module load pyjetty/pyjetty_${build_with_python}
separator "LHAPDF6"
${THISD}/setup_lhapdf6.sh 		--version=6.2.3 	 $@
separator "HEPMC2"
${THISD}/setup_hepmc2_cmake.sh 	--version=2.06.09 	 $@
separator "HEPMC3"
${THISD}/setup_hepmc3.sh 		--version=3.1.1  	 $@
separator "PYTHIA8"
module use ${THISD}/../modules
module load pyjetty/${build_with_python}/HEPMC2/2.06.09 
module load pyjetty/${build_with_python}/HEPMC3/3.1.1
module load pyjetty/${build_with_python}/LHAPDF6/6.2.3
#module load pyjetty/${build_with_python}/fastjet/3.3.2
${THISD}/setup_pythia8.sh 		--version=8235 		 $@
if [ "$(get_opt "install" $@)" == "xyes" ]; then
	separator "redo HEPMC3"
	note "... running with install"
	. ${THISD}/setup_hepmc3.sh 		--version=3.1.1 --re $@
fi
separator "FASTJET"
${THISD}/setup_fastjet.sh 		--version=3.3.2 	 $@
cd ${cdir}
