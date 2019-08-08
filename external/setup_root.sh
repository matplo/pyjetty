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
[ -z ${version} ] && version=6.18.00
note "... version ${version}"
fname=fastjet-${version}
fname=root_v${version}.source
dirsrc=${THISD}/build/${fname}
dirinst=${THISD}/packages/${fname}-${PYJETTY_USER_PYTHON_VERSION}
dirbuild=${dirsrc}-build-${PYJETTY_USER_PYTHON_VERSION}

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
	wget https://root.cern/download/${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${THISD}/build
	tar zxvf ${dirsrc}.tar.gz
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	mkdir -p ${dirbuild}
	if [ -d ${dirsrc} ]; then
		mkdir -p ${dirbuild}
		cd ${dirbuild}
		local _gff=$(which gfortran)
		local _gcc=$(which gcc)
		local _gpp=$(which g++)
		config_opts="-DPYTHON_EXECUTABLE=${PYJETTY_PYTHON_EXECUTABLE} -DPYTHON_INCLUDE_DIR=${PYJETTY_PYTHON_INCLUDE_DIR} -DPYTHON_LIBRARY=${PYJETTY_PYTHON_LIBDIR}"
		config_opts="-Dbuiltin_xrootd=ON -Dmathmore=ON"
		compiler_opts="-DCMAKE_C_COMPILER=${_gcc} -DCMAKE_CXX_COMPILER=${_gpp} -DCMAKE_Fortran_COMPILER=${_gff}"
		echo_info "extra options: ${config_opts} ${compiler_opts}"

		cmake -DCMAKE_BUILD_TYPE\=${BT_build_type} ${compiler_opts} ${config_opts} ${BT_src_dir}

		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && exit 0
		cmake --build . -- -j && cmake -DCMAKE_INSTALL_PREFIX=${dirinst} -P cmake_install.cmake

		cd ${cdir}
	fi
fi

python_path="${dirinst}/lib/ROOT.py"
python_path64="${dirinst}/lib64/ROOT.py"
if [ -f ${python_path64} ] || [ -f ${python_path} ]; then
	${THISD}/../scripts/make_module.sh --dir=${dirinst} --name=ROOT --version=${version}
else
	error "ROOT.py not found - no module generation"
	[ ! -f ${python_path} ] && error "missing ${python_path}"
	[ ! -f ${python_path64} ] && error "missing ${python_path64}"
	exit 0
fi


cd ${cdir}
