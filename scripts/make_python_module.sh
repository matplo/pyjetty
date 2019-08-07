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

# add_path_module "/usr/lib" "USRLIB" ./testmodule
# setenv_module "/usr/lib" "USRLIBSETENV" ./testmodule
# setalias_module $(which python) "userpython" ./testmodule

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

else
	error "no python for ${PYJETTY_USER_PYTHON_VERSION}"
	[ ! -f "${PYJETTY_PYTHON_EXECUTABLE}" ] && error "missing: ${PYJETTY_USER_PYTHON_VERSION}"
	[ ! -f "${PYJETTY_PYTHON_CONFIG_EXECUTABLE}" ] && error "missing: ${PYJETTY_USER_PYTHON_VERSION}-config"
fi
