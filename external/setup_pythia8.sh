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

version=8235
fname=pythia${version}
dirsrc=${SCRIPTPATH}/build/pythia${version}
dirinst=${SCRIPTPATH}/packages/pythia${version}

if [ ! -z ${1} ]; then
	dirinst=${1}
fi

[ ! -d ${SCRIPTPATH}/build ] && mkdir -v ${SCRIPTPATH}/build
[ ! -d ${SCRIPTPATH}/packages ] && mkdir -v ${SCRIPTPATH}/packages

if [ ! -e ${SCRIPTPATH}/build/${fname}.tgz ]; then
	cd ${SCRIPTPATH}/build
	wget http://home.thep.lu.se/~torbjorn/pythia8/${fname}.tgz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fname}.tgz
fi

reconfigure=$([ "x${1}" == "xreconfigure" ] && echo "yes")
# echo "? RECONFIGURE: ${reconfigure}"

function run_configure()
{
	if [ -z ${LHAPDF6_DIR} ]; then
		./configure --prefix=${dirinst} \
			--with-python-include=${python_inc_dir} \
			--with-python-bin=${python_bin_dir} \
			--with-hepmc2=${HEPMC2_DIR} \
			--with-hepmc2-include=${HEPMC2_DIR}/include \
			--with-hepmc2-lib=${HEPMC2_DIR}/lib
	else
		./configure --prefix=${dirinst} \
			--with-python-include=${python_inc_dir} \
			--with-python-bin=${python_bin_dir} \
			--with-hepmc2=${HEPMC2_DIR} \
			--with-hepmc2-include=${HEPMC2_DIR}/include \
			--with-hepmc2-lib=${HEPMC2_DIR}/lib \
			--with-lhapdf6=${LHAPDF6_DIR} \
			--with-lhapdf6-include=${LHAPDF6_DIR}/include \
			--with-lhapdf6-lib=${LHAPDF6_DIR}/lib
	fi
}

if [ ! -d ${dirinst} ] || [ ${reconfigure} ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
	    # echo "unsetting PYTHONPATH"
		[ "x${1}" == "xunset" ] && unset PYTHONPATH	&& echo "unsetting PYTHONPATH"
	    python_inc_dir=$(python3-config --includes | cut -d' ' -f 1 | cut -dI -f 2)
	    python_exec=$(which python3)
	    python_bin_dir=$(dirname ${python_exec})
	    echo "[i] python exec: ${python_exec}"
	    echo "[i] python bin dir: ${python_bin_dir}"
	    echo "[i] python include: ${python_inc_dir}"
	    if [ ! -e ${python_bin_dir}/python ]; then
	    	python_bin_dir=${SCRIPTPATH}/build/pythia-python-bin
	    	mkdir ${python_bin_dir}
	    	ln -s ${python_exec} ${python_bin_dir}/python 
		    echo "[i] fix-up-python bin dir: ${python_bin_dir}"
		fi
    	if [ ! -d ${dirinst} ]; then
			run_configure
			make -j $(n_cores) && make install
		fi
		cd - 2>&1 > /dev/null
	fi
fi

if [ -d ${dirinst} ]; then
	export PYTHIA_DIR=${dirinst}
	export PYTHIA8_ROOT_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi

cd ${cdir}