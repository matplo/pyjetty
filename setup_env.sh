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
STHISDIR=$(thisdir)

if [ -z ${1} ]; then
	source ${STHISDIR}/.config_fjpydev
	echo "[i] FJPYDEV=${FJPYDEV}"
	source ${FJPYDEV}/test/setup_env.sh
else
	echo "[i] using ${1} for fjpydev location"
fi

export PYTHONPATH=${PYTHONPATH}:${STHISDIR}:${STHISDIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${STHISDIR}/lib
