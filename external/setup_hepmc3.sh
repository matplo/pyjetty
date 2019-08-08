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
[ -z ${version} ] && version=3.0.0
note "... version ${version}"
fname=HepMC3-${version}
if [ "x${version}" == "x3.0.0" ]; then
	fname=hepmc${version}
fi
dirsrc=${THISD}/build/${fname}
dirinst=${THISD}/packages/hepmc-${version}-${PYJETTY_USER_PYTHON_VERSION}
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

if [ "x${version}" == "x3.0.0" ]; then
	archsuffix='.tgz'
else
	archsuffix='.tar.gz'
fi

if [ ! -e ${THISD}/build/${fname}${archsuffix} ]; then
	cd ${THISD}/build
	wget http://hepmc.web.cern.ch/hepmc/releases/${fname}${archsuffix}
fi

if [ ! -d ${dirsrc} ]; then
	cd ${THISD}/build
	tar zxvf ${fname}${archsuffix}
fi

redo=$(get_opt "rebuild" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	mkdir -p ${dirbuild}
	if [ -d ${dirsrc} ]; then
		mkdir -p ${dirbuild}
		cd ${dirbuild}
		ROOTIOFLAG=OFF
		[ ! -z ${ROOTSYS} ] && [ -d ${ROOTSYS} ] && ROOTIOFLAG=ON
		echo_info "installing to ${dirinst}"
		cmake \
			-DCMAKE_INSTALL_PREFIX=${dirinst} \
			-DHEPMC3_ENABLE_ROOTIO=${ROOTIOFLAG} \
			-DHEPMC3_BUILD_EXAMPLES=OFF \
			-DHEPMC3_ENABLE_TEST=OFF \
			-DHEPMC3_INSTALL_INTERFACES=ON \
	      	-DCMAKE_MACOSX_RPATH=ON \
	      	-DCMAKE_INSTALL_RPATH=${dirinst}/lib \
	      	-DCMAKE_BUILD_WITH_INSTALL_NAME_DIR=ON \
	      	-DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) \
	      	${dirsrc}
		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && exit 0
		if [ "x${version}" == "x3.0.0" ]; then
			make -j $(n_cores) && make install
		else
			make -j $(n_cores) && make install
			# make test
		fi
		echo_info "link: ln -s ${dirinst}/include/HepMC3 ${dirinst}/include/HepMC"
		ln -s ${dirinst}/include/HepMC3 ${dirinst}/include/HepMC
		echo_info "link ${dirinst}/lib/libHepMC3.dylib"
		_libdir=${dirinst}/lib
		[ -e ${dirinst}/lib64 ] && ln -s ${dirinst}/lib64 ${_libdir} && _libdir=${dirinst}/lib64
		[ -e ${_libdir}/libHepMC3.dylib ] && ln -s ${_libdir}/libHepMC3.dylib ${_libdir}/libHepMC.dylib
		[ -e ${_libdir}/libHepMC3.so ] && ln -s ${_libdir}/libHepMC3.so ${_libdir}/libHepMC.so
		cd ${cdir}
	fi
fi

${THISD}/../scripts/make_module.sh --dir=${dirinst} --name=HEPMC3 --version=${version}

cd ${cdir}
