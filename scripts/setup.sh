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
source ${THISD}/util.sh

if [ ! -f "${HEPPYDIR}/scripts/util.sh" ]; then
    error "this setup relies on HEPPYDIR"
    error "check if modules loaded... "
    exit 0
fi

function make_module()
{
	mkdir -p "${THISD}/../modules/pyjetty"
	modulefiledir=$(abspath_python_expand "${THISD}/../modules/pyjetty")
	module_name="main_${HEPPY_USER_PYTHON_VERSION}"
	modulefile="${modulefiledir}/${module_name}"

	separator "making ${package_name} module ${modulefile}"
	[ -f ${modulefile} ] && warning "removing ${modulefile}" && rm -f ${modulefile}

	heppy_dir=$(abspath_python_expand "${THISD}/..")

	setenv_module ${modulefile} "PYJETTYDIR" ${heppy_dir}
	setenv_module ${modulefile} "PYJETTY_DIR" ${heppy_dir}
	setenv_module ${modulefile} "PYJETTY_ROOT" ${heppy_dir}

	echo "prereq heppy/main_${HEPPY_USER_PYTHON_VERSION}"	>> ${modulefile}

	add_path_module ${modulefile} PYTHONPATH ${heppy_dir}
	setalias_module ${modulefile} heppy_cd "cd ${heppy_dir}"

	dirinst=${1}
	echo_info "installation directory ${dirinst}"
	bin_path="${dirinst}/bin"
	lib_path="${dirinst}/lib"
	lib64_path="${dirinst}/lib64"
	python_path="${dirinst}/lib/python${HEPPY_PYTHON_VERSION}/site-packages"
	python_path64="${dirinst}/lib64/python${HEPPY_PYTHON_VERSION}/site-packages"

	package_name=${2}
	echo_info "package_name ${package_name}"
	setenv_module ${modulefile} ${package_name}DIR ${dirinst}
	setenv_module ${modulefile} ${package_name}_DIR ${dirinst}
	setenv_module ${modulefile} ${package_name}_ROOT ${dirinst}

	setenv_module ${modulefile} ${package_name}_INCLUDE_DIR ${dirinst}/include

	[ $(os_linux) ] && add_path_module ${modulefile} PATH ${bin_path}
	[ $(os_darwin) ] && add_path_module ${modulefile} PATH ${bin_path}

	for sp in ${lib_path} ${lib64_path} ${python_path} ${python_path64}
	do
		[ $(os_linux) ] && add_path_module ${modulefile} LD_LIBRARY_PATH ${sp}
		[ $(os_darwin) ] && add_path_module ${modulefile} DYLD_LIBRARY_PATH ${sp}
	done

	for sp in ${python_path} ${python_path64} ${lib_path} ${lib64_path}
	do
		[ $(os_linux) ] &&  add_path_module ${modulefile} PYTHONPATH ${sp}
		[ $(os_darwin) ] &&  add_path_module ${modulefile} PYTHONPATH ${sp}
	done
}
export -f make_module


separator "${BASH_SOURCE}"

module use ${HEPPYDIR}/modules
module load heppy/main_${HEPPY_USER_PYTHON_VERSION}

${THISD}/../cpptools/scripts/build_cpptools.sh $@

clean=$(get_opt "clean" $@)
cleanall=$(get_opt "cleanall" $@)
if [ -z ${cleanall} ] && [ -z ${clean} ]; then
	make_module $(abspath_python_expand ${THISD}/../cpptools) pyjetty_cpptools
fi

cd ${cpwd}
