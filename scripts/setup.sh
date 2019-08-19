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
source ${THISD}/util.sh
separator "${BASH_SOURCE}"
build_with_python="python"
[ "x$(get_opt "python2" $@)" == "xyes" ] && build_with_python="python2"
[ "x$(get_opt "python3" $@)" == "xyes" ] && build_with_python="python3"
${THISD}/../scripts/make_module.sh --${build_with_python}
module use ${THISD}/../modules
module avail
separator ''
module load pyjetty/pyjetty_${build_with_python}

if [ -z ${PYJETTY_PYTHON_SETUP} ]; then
    error "this setup relies on PYJETTY_PYTHON_SETUP..."
    error "check if modules loaded... module load pyjetty/pyjetty_python ?"
    exit 0
fi

# buildext=$(get_opt "buildext" $@)
# [ "x${buildext}" == "xyes" ] && ${THISD}/../external/setup.sh $@

module load heppy/${build_with_python}/HEPMC2
module load heppy/${build_with_python}/HEPMC3
module load heppy/${build_with_python}/LHAPDF6
module load heppy/${build_with_python}/PYTHIA8
module load heppy/${build_with_python}/FASTJET

build_root=$(get_opt "root" $@)
[ "x${build_root}" == "xyes" ] && module load heppy/${build_with_python}/ROOT

# ( [ ! -d ${THISD}/../cpptools/lib ] || [ "x${redo}" == "xyes" ] ) && ${THISD}/../cpptools/scripts/build_cpptools.sh $@
${THISD}/../cpptools/scripts/build_cpptools.sh $@

clean=$(get_opt "clean" $@)
cleanall=$(get_opt "cleanall" $@)
if [ -z ${cleanall} ] && [ -z ${clean} ]; then
	${THISD}/make_module.sh --make-main-module $@
fi

cd ${cpwd}
