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

version=3.3.2
fjfname=fastjet-${version}
fjdirsrc=${SCRIPTPATH}/build/fastjet-${version}
fjdirinst=${SCRIPTPATH}/packages/fastjet-${version}

if [ ! -z ${1} ]; then
	fjdirinst=${1}
fi

[ ! -d ${SCRIPTPATH}/build ] && mkdir -v ${SCRIPTPATH}/build
[ ! -d ${SCRIPTPATH}/packages ] && mkdir -v ${SCRIPTPATH}/packages

if [ ! -e ${SCRIPTPATH}/build/${fjfname}.tar.gz ]; then
	cd ${SCRIPTPATH}/build
	wget http://fastjet.fr/repo/${fjfname}.tar.gz
fi

if [ ! -d ${fjdirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fjdirsrc}.tar.gz
fi

if [ ! -d ${fjdirinst} ]; then
	if [ -d ${fjdirsrc} ]; then
		cd ${fjdirsrc}
		echo "current dir: $PWD"
		[ "x${1}" == "xunset" ] && unset PYTHONPATH	&& echo "unsetting PYTHONPATH"
	    python_includes=$(python3-config --includes)
	    python_inc_dir=$(python3-config --includes | cut -d' ' -f 1 | cut -dI -f 2)
	    python_exec=$(which python3)
	    python_bin_dir=$(dirname ${python_exec})
	    echo "[i] python exec: ${python_exec}"
	    echo "[i] python includes: ${python_includes}"
	    echo "[i] python include: ${python_inc_dir}"
		if [ "x${CGAL_DIR}" == "x" ]; then
		    ./configure --prefix=${fjdirinst} --enable-allcxxplugins \
		    PYTHON=${python_exec} \
		    PYTHON_INCLUDE="${python_includes}" \
		else
			echo "[i] building using cgal at ${CGAL_DIR}"
		    ./configure --prefix=${fjdirinst} --enable-allcxxplugins \
		    PYTHON=${python_exec} \
		    PYTHON_INCLUDE="${python_includes}" \
		    --enable-cgal --with-cgaldir=${CGAL_DIR} \
		    --enable-pyext
		    # \ LDFLAGS=-Wl,-rpath,${BOOST_DIR}/lib CXXFLAGS=-I${BOOST_DIR}/include CPPFLAGS=-I${BOOST_DIR}/include
		fi
		make -j $(n_cores) && make install
		cd -
	fi
fi

if [ -d ${fjdirinst} ]; then
	export FASTJET_DIR=${fjdirinst}
	export PATH=$PATH:${fjdirinst}/bin
	python_version=$(python3 --version | cut -f 2 -d' ' | cut -f 1-2 -d.)
	export PYTHONPATH=${PYTHONPATH}:${fjdirinst}/lib/python${python_version}/site-packages
fi

cd ${cdir}