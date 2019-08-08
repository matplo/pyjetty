#!/bin/bash

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
export -f thisdir

THISD=$(thisdir)
. ${THISD}/util.sh

function make_python_module()
{
	[ "x$(get_opt "python2" $@)" == "xyes" ] && export PYJETTY_USER_PYTHON_VERSION=python2
	[ "x$(get_opt "python3" $@)" == "xyes" ] && export PYJETTY_USER_PYTHON_VERSION=python3
	[ -z ${PYJETTY_USER_PYTHON_VERSION} ] && export PYJETTY_USER_PYTHON_VERSION=python

	PYJETTY_PYTHON_EXECUTABLE=$(which ${PYJETTY_USER_PYTHON_VERSION})
	PYJETTY_PYTHON_CONFIG_EXECUTABLE=$(which ${PYJETTY_USER_PYTHON_VERSION}-config)
	if [ -f "${PYJETTY_PYTHON_EXECUTABLE}" ] && [ -f "${PYJETTY_PYTHON_CONFIG_EXECUTABLE}" ]; then

		PYJETTY_PYTHON_VERSION=$(${PYJETTY_PYTHON_EXECUTABLE} --version 2>&1 | cut -f 2 -d' ' | cut -f 1-2 -d.)
		PYJETTY_PYTHON_BIN_DIR=$(dirname ${PYJETTY_PYTHON_EXECUTABLE})
		PYJETTY_PYTHON_INCLUDE_DIR=$(${PYJETTY_PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")
		PYJETTY_PYTHON_LIBDIR=$(${PYJETTY_PYTHON_EXECUTABLE} -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
		PYJETTY_PYTHON_NUMPY_INCLUDE_DIR=$(${PYJETTY_PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())")

		PYJETTY_PYTHON_LIBS=$(${PYJETTY_PYTHON_CONFIG_EXECUTABLE} --libs)
		PYJETTY_PYTHON_LIBS_LINK="-L${PYJETTY_PYTHON_LIBDIR} ${PYJETTY_PYTHON_LIBS}"
		PYJETTY_PYTHON_CONFIG_LDFLAGS=$(${PYJETTY_PYTHON_CONFIG_EXECUTABLE} --ldflags)
		PYJETTY_PYTHON_CONFIG_INCLUDES=$(${PYJETTY_PYTHON_CONFIG_EXECUTABLE} --includes)
		PYJETTY_PYTHON_SETUP=TRUE

		modulefiledir=$(abspath ${THISD}/../modules)
		mkdir -p ${modulefiledir}
		modulefile="${modulefiledir}/pyjetty_${PYJETTY_USER_PYTHON_VERSION}"
		separator "making python module ${modulefile}"
		[ -f ${modulefile} ] && warning "removing ${modulefile}" && rm -f ${modulefile}

		setenv_module ${modulefile} PYJETTY_PYTHON_VERSION "${PYJETTY_USER_PYTHON_VERSION}" 
		setenv_module ${modulefile} PYJETTY_PYTHON_EXECUTABLE "${PYJETTY_PYTHON_EXECUTABLE}" 
		setalias_module ${modulefile} pyjetty_show_python "echo ${PYJETTY_USER_PYTHON_VERSION} at ${PYJETTY_PYTHON_EXECUTABLE}" 
		setenv_module ${modulefile} PYJETTY_USER_PYTHON_VERSION "${PYJETTY_USER_PYTHON_VERSION}"

		setenv_module ${modulefile} PYJETTY_PYTHON_VERSION ${PYJETTY_PYTHON_VERSION}
		setenv_module ${modulefile} PYJETTY_PYTHON_BIN_DIR ${PYJETTY_PYTHON_BIN_DIR} 
		setenv_module ${modulefile} PYJETTY_PYTHON_INCLUDE_DIR ${PYJETTY_PYTHON_INCLUDE_DIR}
		setenv_module ${modulefile} PYJETTY_PYTHON_LIBDIR ${PYJETTY_PYTHON_LIBDIR} 
		setenv_module ${modulefile} PYJETTY_PYTHON_NUMPY_INCLUDE_DIR ${PYJETTY_PYTHON_NUMPY_INCLUDE_DIR} 

		setenv_module ${modulefile} PYJETTY_PYTHON_LIBS ${PYJETTY_PYTHON_LIBS} 
		setenv_module ${modulefile} PYJETTY_PYTHON_LIBS_LINK ${PYJETTY_PYTHON_LIBS_LINK} 
		setenv_module ${modulefile} PYJETTY_PYTHON_CONFIG_LDFLAGS ${PYJETTY_PYTHON_CONFIG_LDFLAGS} 
		setenv_module ${modulefile} PYJETTY_PYTHON_CONFIG_INCLUDES ${PYJETTY_PYTHON_CONFIG_INCLUDES} 
		setenv_module ${modulefile} PYJETTY_PYTHON_SETUP ${PYJETTY_PYTHON_SETUP} 

		setenv_module ${modulefile} PYJETTY_PYTHON_MODULE_LOADED "pyjetty_${PYJETTY_USER_PYTHON_VERSION}"

	else
		error "no python for ${PYJETTY_USER_PYTHON_VERSION}"
		[ ! -f "${PYJETTY_PYTHON_EXECUTABLE}" ] && error "missing: ${PYJETTY_USER_PYTHON_VERSION}"
		[ ! -f "${PYJETTY_PYTHON_CONFIG_EXECUTABLE}" ] && error "missing: ${PYJETTY_USER_PYTHON_VERSION}-config"
	fi
}
export -f make_python_module

function make_module_package()
{
	dirinst=${1}
	package_name=$(basename ${dirinst})
	python_path="${dirinst}/lib/python${PYJETTY_PYTHON_VERSION}/site-packages/fastjet.py"
	python_path64="${dirinst}/lib64/python${PYJETTY_PYTHON_VERSION}/site-packages/fastjet.py"

	if [ -f ${python_path64} ] || [ -f ${python_path} ]; then
		modulefiledir=$(abspath ${THISD}/../modules)
		mkdir -p ${modulefiledir}
		modulefile="${modulefiledir}/pyjetty_${package_name}"
		separator "making ${package_name} module ${modulefile}"
		[ -f ${modulefile} ] && warning "removing ${modulefile}" && rm -f ${modulefile}

		bin_path="${dirinst}/bin"
		lib_path="${dirinst}/lib"
		lib64_path="${dirinst}/lib64"
		python_path="${dirinst}/lib/python${PYJETTY_PYTHON_VERSION}/site-packages"
		python_path64="${dirinst}/lib64/python${PYJETTY_PYTHON_VERSION}/site-packages"

		setenv_module ${modulefile} ${package_name}DIR ${dirinst}
		setenv_module ${modulefile} ${package_name}_DIR ${dirinst}
		setenv_module ${modulefile} ${package_name}_ROOT ${dirinst}

		[ $(os_linux) ] && add_path_module ${modulefile} PATH ${bin_path}
		[ $(os_darwin) ] && add_path_module ${modulefile} PATH ${bin_path}

		for sp in ${lib_path} ${lib64_path} ${python_path} ${python_path64}
		do
			[ $(os_linux) ] && add_path_module ${modulefile} LD_LIBRARY_PATH ${sp}
			[ $(os_darwin) ] && add_path_module ${modulefile} DYLD_LIBRARY_PATH ${sp}
		done

		for sp in ${python_path} ${python_path64}
		do
			[ $(os_linux) ] && add_path_module ${modulefile} PYTHONPATH ${sp}
		done

		if [ ! -z ${PYJETTY_PYTHON_MODULE_LOADED} ]; then
			echo "module prereq ${python_module_name}" >> ${modulefile}
		fi

	else
		error "fastjet.py not found - no module generation"
		[ ! -f ${python_path} ] && error "missing ${python_path}"
		[ ! -f ${python_path64} ] && error "missing ${python_path64}"
	fi
}
export -f make_module_package

if [ "x$(get_opt "python" $@)" == "xyes" ] || [ "x$(get_opt "python2" $@)" == "xyes" ] || [ "x$(get_opt "python2" $@)" == "xyes" ]; then
	separator "make_modules.sh :: python module"
	make_python_module $@
	separator "make_modules.sh - done"
fi

packagedir=$(get_opt "packagedir" $@)
if [ -d "${packagedir}" ]; then
	separator "make_modules.sh :: package module"
	make_module_package ${packagedir}
	separator "make_modules.sh - done"
fi
