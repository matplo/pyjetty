#!/bin/bash

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

version=2.06.09
fname=HepMC-${version}
dirsrc=${SCRIPTPATH}/build/HepMC-${version}
dirinst=${SCRIPTPATH}/packages/hepmc-${version}

BT_remote_file=http://lcgapp.cern.ch/project/simu/HepMC/download/HepMC-${BT_version}.tar.gz

if [ ! -z ${1} ]; then
	dirinst=${1}
fi

[ ! -d ${SCRIPTPATH}/build ] && mkdir -v ${SCRIPTPATH}/build
[ ! -d ${SCRIPTPATH}/packages ] && mkdir -v ${SCRIPTPATH}/packages

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
		[ "x${1}" == "xunset" ] && unset PYTHONPATH	&& echo "unsetting PYTHONPATH"
	    python_inc_dir=$(python3-config --includes | cut -d' ' -f 1 | cut -dI -f 2)
	    python_exec=$(which python3)
	    python_bin_dir=$(dirname ${python_exec})
	    echo "[i] python exec: ${python_exec}"
	    echo "[i] python include: ${python_inc_dir}"
	    # this is a nasty trick to force python3 bindings
	    python_bin_dir="$PWD/tmppy"
	    mkdir -p ${python_bin_dir}
	    ln -s ${python_exec} ${python_bin_dir}/python
	    echo "python bin dir: ${python_bin_dir}"
		./configure --prefix=${dirinst} --with-momentum=GEV --with-length=CM
		make -j $(n_cores) && make install
		cd ${savepwd}
	fi
fi

if [ -d ${dirinst} ]; then
	export HEPMC2_DIR=${dirinst}
	export HEPMC_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi
