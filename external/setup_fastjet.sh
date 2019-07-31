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
[ -z ${version} ] && version=3.3.2
note "... version ${version}"
fname=fastjet-${version}
dirsrc=${SCRIPTPATH}/build/fastjet-${version}
dirinst=${SCRIPTPATH}/packages/fastjet-${version}

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
	wget http://fastjet.fr/repo/${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${dirsrc}.tar.gz
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
		[ "$(get_opt "unset" $@)" == "xyes" ] && unset PYTHONPATH && warning "unsetting PYTHONPATH"
	    python_includes=$(python3-config --includes)
	    python_inc_dir=$(python3-config --includes | cut -d' ' -f 1 | cut -dI -f 2)
	    python_exec=$(which python3)
	    python_bin_dir=$(dirname ${python_exec})
	    echo_info "python exec: ${python_exec}"
	    echo_info "python includes: ${python_includes}"
	    echo_info "python include: ${python_inc_dir}"
	    if [ "x${CGAL_DIR}" == "x" ] || [ ! -d ${CGAL_DIR} ]; then
    		no_cgal=$(get_opt "no-cgal" $@)
	    	if [ "x$(os_darwin)" == "xyes" ] && [ -z ${no_cgal} ]; then
	    		[ -e /usr/local/lib/libCGAL.dylib ] && [ -e /usr/local/include/CGAL/config.h ] && export CGAL_DIR=/usr/local
	    		note "trying CGAL in ${CGAL_DIR} - override with --no-cgal"
	    	fi
	    fi
		if [ "x${CGAL_DIR}" == "x" ] || [ ! -d ${CGAL_DIR} ]; then
			note "building w/o cgal"
		    ./configure --prefix=${dirinst} --enable-allcxxplugins \
		    PYTHON=${python_exec} \
		    PYTHON_INCLUDE="${python_includes}" \
		    --enable-pyext
		else
			note "building using cgal at ${CGAL_DIR}"
		    ./configure --prefix=${dirinst} --enable-allcxxplugins \
		    PYTHON=${python_exec} \
		    PYTHON_INCLUDE="${python_includes}" \
		    --enable-cgal --with-cgaldir=${CGAL_DIR} \
		    --enable-pyext
		    # \ LDFLAGS=-Wl,-rpath,${BOOST_DIR}/lib CXXFLAGS=-I${BOOST_DIR}/include CPPFLAGS=-I${BOOST_DIR}/include
		fi
		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && return 0
		make -j $(n_cores) && make install
		cd ${cdir}
	fi
fi

if [ -d ${dirinst} ]; then
	export FASTJET_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	python_version=$(python3 --version | cut -f 2 -d' ' | cut -f 1-2 -d.)
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib/python${python_version}/site-packages
	[ -d ${dirinst}/lib64 ] && export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib64/python${python_version}/site-packages
fi

cd ${cdir}
