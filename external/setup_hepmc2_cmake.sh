#!/bin/bash

cdir=$(pwd)

function os_linux()
{
	_system=$(uname -a | cut -f 1 -d " ")
	if [ $_system == "Linux" ]; then
		echo "yes"
	else
		echo
	fi
}

function os_darwin()
{
	_system=$(uname -a | cut -f 1 -d " ")
	if [ $_system == "Darwin" ]; then
		echo "yes"
	else
		echo
	fi
}

function n_cores()
{
	local _ncores="1"
	[ $(os_darwin) ] && local _ncores=$(system_profiler SPHardwareDataType | grep "Number of Cores" | cut -f 2 -d ":" | sed 's| ||')
	[ $(os_linux) ] && local _ncores=$(lscpu | grep "CPU(s):" | head -n 1 | cut -f 2 -d ":" | sed 's| ||g')
	#[ ${_ncores} -gt "1" ] && retval=$(_ncores-1)
	echo ${_ncores}
}

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

savepwd=${PWD}

version=2.06.10
fname=HepMC-${version}
dirsrc=${SCRIPTPATH}/build/HepMC-${version}
dirinst=${SCRIPTPATH}/packages/hepmc-${version}

[ ! -d ${SCRIPTPATH}/build ] && mkdir -v ${SCRIPTPATH}/build
[ ! -d ${SCRIPTPATH}/packages ] && mkdir -v ${SCRIPTPATH}/packages

if [ ! -z ${1} ]; then
	dirinst=${1}
fi

if [ ! -e ${SCRIPTPATH}/build/${fname}.tar.gz ]; then
	cd ${SCRIPTPATH}/build
	wget http://lcgapp.cern.ch/project/simu/HepMC/download/${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fname}.tar.gz
fi

if [ ! -d ${dirinst} ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
		[ "x${version}" == "x2.06.09" ] && patch CMakeLists.txt -i ${SCRIPTPATH}/patches/HepMC-2.06.09-CMakeLists.txt.patch
		mkdir ${SCRIPTPATH}/build/build_dir_${fname}
		cd ${SCRIPTPATH}/build/build_dir_${fname}
		cmake -Dmomentum:STRING=GEV -Dlength:STRING=CM \
				-DCMAKE_INSTALL_PREFIX=${dirinst} \
		     	-DCMAKE_BUILD_TYPE=Release \
		      	-Dbuild_docs:BOOL=OFF \
		      	-DCMAKE_MACOSX_RPATH=ON \
		      	-DCMAKE_INSTALL_RPATH=${dirinst}/lib \
		      	-DCMAKE_BUILD_WITH_INSTALL_NAME_DIR=ON \
			    ${dirsrc}
		make && make install
		make test
		cd ${savepwd}
	fi
fi


if [ -d ${dirinst} ]; then
	export HEPMC2_DIR=${dirinst}
	export HEPMC_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi

cd ${cdir}