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
[ "x$PYJETTY_PYTHON_SETUP" != "xTRUE" ] && error "missing: PYJETTY_PYTHON_SETUP" && exit 1
warning "user python version: ${PYJETTY_USER_PYTHON_VERSION}"
version=$(get_opt "version" $@)
[ -z ${version} ] && version=3.3.2
note "... version ${version}"
fname=fastjet-${version}
dirsrc=${THISD}/build/fastjet-${version}
warning "user python version: ${PYJETTY_USER_PYTHON_VERSION}"
dirinst=${THISD}/packages/fastjet-${version}-${PYJETTY_USER_PYTHON_VERSION}

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
	wget http://fastjet.fr/repo/${fname}.tar.gz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${THISD}/build
	tar zxvf ${dirsrc}.tar.gz
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
	    if [ "x${CGAL_DIR}" == "x" ] || [ ! -d ${CGAL_DIR} ]; then
    		no_cgal=$(get_opt "no-cgal" $@)
	    	if [ "x$(os_darwin)" == "xyes" ] && [ -z ${no_cgal} ]; then
	    		[ -e /usr/local/lib/libCGAL.dylib ] && [ -e /usr/local/include/CGAL/config.h ] && export CGAL_DIR=/usr/local
	    		note "trying CGAL in ${CGAL_DIR} - override with --no-cgal"
	    	fi
	    fi
		if [ "x${CGAL_DIR}" == "x" ] || [ ! -d ${CGAL_DIR} ]; then
			note "building w/o cgal"
		    PYTHON=${PYJETTY_PYTHON_EXECUTABLE} \
		    PYTHON_INCLUDE="${PYJETTY_PYTHON_INCLUDES}" \
		    ./configure --prefix=${dirinst} --enable-allcxxplugins \
		    --enable-pyext
		else
			note "building using cgal at ${CGAL_DIR}"
		    PYTHON=${PYJETTY_PYTHON_EXECUTABLE} \
		    PYTHON_INCLUDE="${PYJETTY_PYTHON_INCLUDES}" \
		    ./configure --prefix=${dirinst} --enable-allcxxplugins \
		    --enable-cgal --with-cgaldir=${CGAL_DIR} \
		    --enable-pyext
		    # \ LDFLAGS=-Wl,-rpath,${BOOST_DIR}/lib CXXFLAGS=-I${BOOST_DIR}/include CPPFLAGS=-I${BOOST_DIR}/include
		fi
		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && exit 0
		make -j $(n_cores) && make install
		cd ${cdir}
	fi
fi

python_path="${dirinst}/lib/python${PYJETTY_PYTHON_VERSION}/site-packages/fastjet.py"
python_path64="${dirinst}/lib64/python${PYJETTY_PYTHON_VERSION}/site-packages/fastjet.py"

if [ -f ${python_path64} ] || [ -f ${python_path} ]; then
	# export FASTJET_DIR=${dirinst}
	modulefiledir=$(abspath ${THISD}/../modules)
	mkdir -p ${modulefiledir}
	modulefile="${modulefiledir}/pyjetty_fastjet_${PYJETTY_USER_PYTHON_VERSION}"
	separator "making python module ${modulefile}"
	[ -f ${modulefile} ] && warning "removing ${modulefile}" && rm -f ${modulefile}

	bin_path="${dirinst}/bin"
	lib_path="${dirinst}/lib"
	lib64_path="${dirinst}/lib64"
	python_path="${dirinst}/lib/python${PYJETTY_PYTHON_VERSION}/site-packages"
	python_path64="${dirinst}/lib64/python${PYJETTY_PYTHON_VERSION}/site-packages"

	setenv_module ${modulefile} FASTJETDIR ${dirinst}
	setenv_module ${modulefile} FASTJET_DIR ${dirinst}
	setenv_module ${modulefile} FASTJET_ROOT ${dirinst}

	for sp in ${lib_path} ${lib64_path} ${python_path} ${python_path64}
	do
		[ $(os_linux) ] && add_path_module ${modulefile} LD_LIBRARY_PATH ${sp}
		[ $(os_darwin) ] && add_path_module ${modulefile} DYLD_LIBRARY_PATH ${sp}
	done

	for sp in ${python_path} ${python_path64}
	do
		[ $(os_linux) ] && add_path_module ${modulefile} PYTHONPATH ${sp}
	done
else
	error "fastjet.py not found - no module generation"
	[ ! -f ${python_path} ] && error "missing ${python_path}"
	[ ! -f ${python_path64} ] && error "missing ${python_path64}"
fi

cd ${cdir}
