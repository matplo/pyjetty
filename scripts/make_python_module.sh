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

function add_path_module()
{
	path=${1}
	what=${2}
	modulefile=${3}
	if [ ! -z "${path}" ] && [ -d ${path} ]; then
		if [ ! -f "${modulefile}" ]; then
			touch ${modulefile}
			if [ ! -f "${modulefile}" ]; then
				error "add_path_module:: unable to open file ${modulefile}"
			else
				echo "#%Module" >> ${modulefile}
				echo "proc ModulesHelp { } {" >> ${modulefile}
				echo "    global version" >> ${modulefile}
				echo "    puts stderr \"   Setup pyjetty ${version}\"}" >> ${modulefile}
				echo "set version ${modulefile}" >> ${modulefile}
			fi
		fi
		if [ ! -f "${modulefile}" ]; then
			error "add_path_module:: unable to open file ${modulefile}"
		else
			echo_info "adding ${what} to ${path}"
			echo "prepend-path ${what} ${path}" >> ${modulefile}
		fi
	else
		error "add_path_module:: ignoring ${path} for ${what}"
	fi
}
export -f add_path_module

function setenv_module()
{
	path=${1}
	what=${2}
	modulefile=${3}
	if [ ! -z "${path}" ]; then
		if [ ! -f "${modulefile}" ]; then
			touch ${modulefile}
			if [ ! -f "${modulefile}" ]; then
				error "setenv_module:: unable to open file ${modulefile}"
			else
				echo "#%Module" >> ${modulefile}
				echo "proc ModulesHelp { } {" >> ${modulefile}
				echo "    global version" >> ${modulefile}
				echo "    puts stderr \"   Setup pyjetty ${version}\"}" >> ${modulefile}
				echo "set version ${modulefile}" >> ${modulefile}
			fi
		fi
		if [ ! -f "${modulefile}" ]; then
			error "setenv_module:: unable to open file ${modulefile}"
		else
			echo_info "setenv ${what} ${path}"
			echo "setenv ${what} ${path}" >> ${modulefile}
		fi
	else
		error "setenv_module:: ignoring ${path} for ${what}"
	fi
}
export -f setenv_module

function setalias_module()
{
	path=${1}
	what=${2}
	modulefile=${3}
	if [ ! -f "${modulefile}" ]; then
		touch ${modulefile}
		if [ ! -f "${modulefile}" ]; then
			error "setalias_module:: unable to open file ${modulefile}"
		else
			echo "#%Module" >> ${modulefile}
			echo "proc ModulesHelp { } {" >> ${modulefile}
			echo "    global version" >> ${modulefile}
			echo "    puts stderr \"   Setup pyjetty ${version}\"}" >> ${modulefile}
			echo "set version ${modulefile}" >> ${modulefile}
		fi
	fi
	if [ ! -f "${modulefile}" ]; then
		error "setalias_module:: unable to open file ${modulefile}"
	else
		echo_info "set-alias ${what} ${path}"
		echo "set-alias ${what} \"${path}\"" >> ${modulefile}
	fi
}
export -f setalias_module

# add_path_module "/usr/lib" "USRLIB" ./testmodule
# setenv_module "/usr/lib" "USRLIBSETENV" ./testmodule
# setalias_module $(which python) "userpython" ./testmodule

[ "x$(get_opt "python2" $@)" == "xyes" ] && export PYJETTY_USER_PYTHON_VERSION=python2
[ "x$(get_opt "python3" $@)" == "xyes" ] && export PYJETTY_USER_PYTHON_VERSION=python3
[ -z ${PYJETTY_USER_PYTHON_VERSION} ] && export PYJETTY_USER_PYTHON_VERSION=python

PYJETTY_PYTHON_EXECUTABLE=$(which ${PYJETTY_USER_PYTHON_VERSION})

if [ -f "${PYJETTY_PYTHON_EXECUTABLE}" ]; then
	modulefiledir=$(abspath ${THISD}/../modules)
	mkdir -p ${modulefiledir}
	modulefile="${modulefiledir}/pyjetty_${PYJETTY_USER_PYTHON_VERSION}"
	rm ${modulefile}
	separator "making python module ${modulefile}"
	setenv_module "$PYJETTY_USER_PYTHON_VERSION" PYJETTY_PYTHON_VERSION ${modulefile}
	setenv_module ${PYJETTY_PYTHON_EXECUTABLE} PYJETTY_PYTHON_EXECUTABLE ${modulefile}
	setenv_module ${PYJETTY_PYTHON_EXECUTABLE} PYJETTY_PYTHON_EXECUTABLE ${modulefile}
	setalias_module "echo ${PYJETTY_USER_PYTHON_VERSION} at ${PYJETTY_PYTHON_EXECUTABLE}" pyjetty_show_python ${modulefile}
else
	error "no python for ${PYJETTY_USER_PYTHON_VERSION}"
fi
