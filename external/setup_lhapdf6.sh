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

npdfs_link="http://lhapdfsets.web.cern.ch/lhapdfsets/current/EPPS16nlo_CT14nlo_Pb208.tar.gz"

version=$(get_opt "version" $@)
[ -z ${version} ] && version=6.2.3
note "... version ${version}"
fname=LHAPDF-${version}
dirsrc=${SCRIPTPATH}/build/LHAPDF-${version}
dirinst=${SCRIPTPATH}/packages/LHAPDF-${version}

function grace_return()
{
	cd ${cdir}
}
prefix=$(get_opt "prefix" $@)
[ ! -z ${prefix} ] && dirinst=${prefix}
clean=$(get_opt "clean" $@)
if [ "x${clean}" == "xyes" ]; then
	warning "cleaning..."
	echo_info "removing ${dirsrc}"
	rm -rf ${dirsrc}
	echo_info "removing ${dirinst}"
	rm -rf ${dirinst}
	grace_return && return 0
fi
uninstall=$(get_opt "uninstall" $@)
if [ "x${uninstall}" == "xyes" ]; then
	echo_info "uninstall..."
	rm -rf ${dirinst}
	grace_return && return 0
fi
installed=$(get_opt "installed" $@)
if [ "x${installed}" == "xyes" ]; then
	[ -d ${dirinst} ] && echo_info "${dirinst} exists"
	[ ! -d ${dirinst} ] && error "${dirinst} does NOT exists"
	grace_return && return 0
fi

[ ! -d ${SCRIPTPATH}/build ] && mkdir -v ${SCRIPTPATH}/build
[ ! -d ${SCRIPTPATH}/packages ] && mkdir -v ${SCRIPTPATH}/packages

if [ ! -e ${SCRIPTPATH}/build/${fname}.tar.gz ]; then
	cd ${SCRIPTPATH}/build
	wget https://lhapdf.hepforge.org/downloads/?f=${fname}.tar.gz -O ${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fname}.tar.gz
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
		./configure --prefix=${dirinst}
		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && return 0
		make -j $(n_cores) && make install
		cd -
	fi
fi

if [ -d ${dirinst} ]; then
	export LHAPDF6_DIR=${dirinst}
	export LHAPDF6_ROOT_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	python_version=$(python3 --version | cut -f 2 -d' ' | cut -f 1-2 -d.)
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib/python${python_version}/site-packages
	export LHAPATH=${dirinst}/share/LHAPDF
fi

cd ${cdir}
