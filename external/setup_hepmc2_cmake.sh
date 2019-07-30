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

version=$(get_opt "version" $@)
[ -z ${version} ] && version=2.06.09
warning "... version ${version}"
fname=HepMC-${version}
dirsrc=${SCRIPTPATH}/build/HepMC-${version}
dirinst=${SCRIPTPATH}/packages/hepmc-${version}

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
	wget http://lcgapp.cern.ch/project/simu/HepMC/download/${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fname}.tar.gz
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
		[ "x${version}" == "x2.06.09" ] && patch -N CMakeLists.txt -i ${SCRIPTPATH}/patches/HepMC-2.06.09-CMakeLists.txt.patch
		mkdir ${SCRIPTPATH}/build/build_dir_${fname}
		cd ${SCRIPTPATH}/build/build_dir_${fname}
		cmake -Dmomentum:STRING=GEV -Dlength:STRING=CM \
				-DCMAKE_INSTALL_PREFIX=${dirinst} \
		     	-DCMAKE_BUILD_TYPE=Release \
		      	-Dbuild_docs:BOOL=OFF \
		      	-DCMAKE_MACOSX_RPATH=ON \
		      	-DCMAKE_INSTALL_RPATH=${dirinst}/lib \
		      	-DCMAKE_BUILD_WITH_INSTALL_NAME_DIR=ON \
		      	-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) \
			    ${dirsrc}
		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && return 0
		make && make install
		make test
		cd ${cdir}
	fi
fi


if [ -d ${dirinst} ]; then
	export HEPMC2_DIR=${dirinst}
	export HEPMC2_ROOT_DIR=${dirinst}
	# export HEPMC_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi

cd ${cdir}