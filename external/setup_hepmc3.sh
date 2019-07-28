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

version=3.1.1
fname=HepMC3-${version}
dirsrc=${SCRIPTPATH}/build/HepMC3-${version}
dirinst=${SCRIPTPATH}/packages/hepmc-${version}

[ ! -d ${SCRIPTPATH}/build ] && mkdir -v ${SCRIPTPATH}/build
[ ! -d ${SCRIPTPATH}/packages ] && mkdir -v ${SCRIPTPATH}/packages

if [ ! -z ${1} ]; then
	dirinst=${1}
fi

if [ ! -e ${SCRIPTPATH}/build/${fname}.tar.gz ]; then
	cd ${SCRIPTPATH}/build
	wget http://hepmc.web.cern.ch/hepmc/releases/${fname}.tar.gz
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
	    echo "python exec: ${python_exec}"
	    echo "python include: ${python_inc_dir}"
	    # this is a nasty trick to force python3 bindings
	    if [ ! -e ${python_bin_dir}/python ]; then
	    	python_bin_dir=${SCRIPTPATH}/build/pythia-python-bin
	    	mkdir ${python_bin_dir}
	    	ln -s ${python_exec} ${python_bin_dir}/python 
		    echo "[i] fix-up-python bin dir: ${python_bin_dir}"
		fi
		mkdir hepmc3-build
		cd hepmc3-build
		cmake -S${dirsrc} -DHEPMC3_ENABLE_ROOTIO=OFF -DCMAKE_INSTALL_PREFIX=${dirinst} \
			-DHEPMC3_BUILD_EXAMPLES=ON -DHEPMC3_ENABLE_TEST=ON \
			-DHEPMC3_INSTALL_INTERFACES=ON \
	      	-DCMAKE_MACOSX_RPATH=ON \
	      	-DCMAKE_INSTALL_RPATH=${dirinst}/lib \
	      	-DCMAKE_BUILD_WITH_INSTALL_NAME_DIR=ON 
		make -j $(n_cores) && make install && make test
		ln -s ${dirinst}/include/HepMC3 ${dirinst}/include/HepMC
		[ -e ${dirinst}/lib/libHepMC3.dylib ] && ln -s ${dirinst}/lib/libHepMC3.dylib ${dirinst}/lib/libHepMC.dylib
		[ -e ${dirinst}/lib/libHepMC3.so ] && ln -s ${dirinst}/lib/libHepMC3.so ${dirinst}/lib/libHepMC.so
		cd ${savepwd}
	fi
fi

if [ -d ${dirinst} ]; then
	export HEPMC3_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi

cd ${cdir}