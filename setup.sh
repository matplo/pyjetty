#!/bin/bash

cpwd=$(pwd)

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
STHISDIR=$(thisdir)
source ${STHISDIR}/external/util.sh
separator "pyjetty $(abspath ${BASH_SOURCE})"
if [ ! -z ${PYJETTY_SET} ]; then
	warning "PYJETTY_SET is ${PYJETTY_SET}"
fi

if [ -z ${PYJETTY_SETUP_EXTERNAL} ]; then
	if [ -e ${PWD}/.pyjetty_config_external ]; then
		source ${PWD}/.pyjetty_config_external
	else
		if [ -e ${STHISDIR}/.pyjetty_config_external ]; then
			source ${STHISDIR}/.pyjetty_config_external
			[ -e ${PYJETTY_SETUP_EXTERNAL} ] && echo_info "[i] PYJETTY_SETUP_EXTERNAL=${PYJETTY_SETUP_EXTERNAL}" && source ${PYJETTY_SETUP_EXTERNAL}
		else
			echo_info "PYJETTY_SETUP_EXTERNAL=${STHISDIR}/external/setup.sh" | tee ${STHISDIR}/.pyjetty_config_external
			if [ -e ${STHISDIR}/.pyjetty_config_external ]; then
				source ${STHISDIR}/.pyjetty_config_external
			fi
		fi
	fi
	[ -e ${PYJETTY_SETUP_EXTERNAL} ] && echo_info "[i] PYJETTY_SETUP_EXTERNAL=${PYJETTY_SETUP_EXTERNAL}" && source ${PYJETTY_SETUP_EXTERNAL}
fi

redo=$(get_opt "re" $@)
( [ ! -d ${STHISDIR}/cpptools/lib ] || [ "x${redo}" == "xyes" ] ) && ${STHISDIR}/cpptools/scripts/build_cpptools.sh $@

export PYJETTY_SET=TRUE
if [ -z ${PYJETTY_SET} ]; then
	export PYTHONPATH=${PYTHONPATH}:${STHISDIR}:${STHISDIR}/cpptools/lib
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${STHISDIR}/cpptools/lib
	export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${STHISDIR}/cpptools/lib
fi

cd ${cpwd}