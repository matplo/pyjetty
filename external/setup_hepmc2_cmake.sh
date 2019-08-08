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
. ${THISD}/../scripts/util.sh
separator "${BASH_SOURCE}"
if [ -z "${PYJETTY_USER_PYTHON_VERSION}" ]; then
	warning "trying to load pyjetty_python"
	pyjetty_python_module_name=$(module avail -t | grep pyjetty_python | head -n 1 | grep pyjetty_python)
	if [ ! -z ${pyjetty_python_module_name} ]; then
		warning "... found ${pyjetty_python_module_name}"	
		module load ${pyjetty_python_module_name}
	else
		warning "... no suitable module found"
	fi
fi
[ -z "${PYJETTY_USER_PYTHON_VERSION}" ] && error "missing: PYJETTY_USER_PYTHON_VERSION" && exit 1
warning "using pyjetty python version: ${PYJETTY_USER_PYTHON_VERSION}"
version=$(get_opt "version" $@)
[ -z ${version} ] && version=2.06.09
note "... version ${version}"
fname=HepMC-${version}
dirsrc=${THISD}/build/HepMC-${version}
dirinst=${THISD}/packages/hepmc-${version}-${PYJETTY_USER_PYTHON_VERSION}

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
	grace_return && exit 0
fi
uninstall=$(get_opt "uninstall" $@)
if [ "x${uninstall}" == "xyes" ]; then
	echo_info "uninstall..."
	rm -rf ${dirinst}
	grace_return && exit 0
fi
installed=$(get_opt "installed" $@)
if [ "x${installed}" == "xyes" ]; then
	[ -d ${dirinst} ] && echo_info "${dirinst} exists"
	[ ! -d ${dirinst} ] && error "${dirinst} does NOT exists"
	grace_return && exit 0
fi

[ ! -d ${THISD}/build ] && mkdir -v ${THISD}/build
[ ! -d ${THISD}/packages ] && mkdir -v ${THISD}/packages

if [ ! -e ${THISD}/build/${fname}.tar.gz ]; then
	cd ${THISD}/build
	wget http://lcgapp.cern.ch/project/simu/HepMC/download/${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${THISD}/build
	tar zxvf ${fname}.tar.gz
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
		[ "x${version}" == "x2.06.09" ] && patch -N CMakeLists.txt -i ${THISD}/patches/HepMC-2.06.09-CMakeLists.txt.patch
		mkdir ${THISD}/build/build_dir_${fname}
		cd ${THISD}/build/build_dir_${fname}
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
		[ "x${configure_only}" == "xyes" ] && grace_return && exit 0
		make && make install
		make test
		cd ${cdir}
	fi
fi

${THISD}/../scripts/make_module.sh --dir=${dirinst} --name=HEPMC2 --version=${version}

cd ${cdir}
