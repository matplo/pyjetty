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
. ${THISD}/util.sh

mdir=$(abspath_python_expand "${THISD}/..")
separator "rm -rf ${mdir}/modules"
rm -rf ${mdir}/modules
#separator "rm -rf ${mdir}/external/build"
#rm -rf ${mdir}/external/build
separator "cleaning ${mdir}/external/build - leaving *.gz files"
find ${mdir}/external/build -type f ! -name '*.gz' -delete
separator "rm -rf ${mdir}/external/packages"
rm -rf ${mdir}/external/packages

cd ${cdir}
