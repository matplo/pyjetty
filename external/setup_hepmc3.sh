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
[ -z ${version} ] && version=3.0.0
warning "... version ${version}"
fname=HepMC3-${version}
if [ "x${version}" == "x3.0.0" ]; then
	fname=hepmc${version}
fi
dirsrc=${SCRIPTPATH}/build/${fname}
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

if [ "x${version}" == "x3.0.0" ]; then
	archsuffix='.tgz'
else
	archsuffix='.tar.gz'
fi

if [ ! -e ${SCRIPTPATH}/build/${fname}${archsuffix} ]; then
	cd ${SCRIPTPATH}/build
	wget http://hepmc.web.cern.ch/hepmc/releases/${fname}${archsuffix}
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fname}${archsuffix}
fi

redo=$(get_opt "re" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
		[ "x${1}" == "xunset" ] && unset PYTHONPATH	&& echo "unsetting PYTHONPATH"
	    python_inc_dir=$(python3-config --includes | cut -d' ' -f 1 | cut -dI -f 2)
	    python_exec=$(which python3)
	    python_bin_dir=$(dirname ${python_exec})
	    echo_info "python exec: ${python_exec}"
	    echo_info "python include: ${python_inc_dir}"
	    # this is a nasty trick to force python3 bindings
	    if [ ! -e ${python_bin_dir}/python ]; then
	    	python_bin_dir=${SCRIPTPATH}/build/pythia-python-bin
	    	mkdir ${python_bin_dir}
	    	ln -s ${python_exec} ${python_bin_dir}/python 
		    warning "fix-up-python bin dir: ${python_bin_dir}"
		fi
		mkdir hepmc3-build
		cd hepmc3-build
		echo_info "installing to ${dirinst}"
		cmake -S${dirsrc} \
			-DCMAKE_INSTALL_PREFIX=${dirinst} \
			-DHEPMC3_ENABLE_ROOTIO=OFF \
			-DHEPMC3_BUILD_EXAMPLES=ON \
			-DHEPMC3_ENABLE_TEST=ON \
			-DHEPMC3_INSTALL_INTERFACES=ON \
	      	-DCMAKE_MACOSX_RPATH=ON \
	      	-DCMAKE_INSTALL_RPATH=${dirinst}/lib \
	      	-DCMAKE_BUILD_WITH_INSTALL_NAME_DIR=ON	    
		configure_only=$(get_opt "configure-only" $@)
		[ "x${configure_only}" == "xyes" ] && grace_return && return 0
		if [ "x${version}" == "x3.0.0" ]; then
			make -j $(n_cores) && make install
		else
			make -j $(n_cores) && make install 
			# make test
		fi
		echo_info "link: ln -s ${dirinst}/include/HepMC3 ${dirinst}/include/HepMC"
		ln -s ${dirinst}/include/HepMC3 ${dirinst}/include/HepMC
		echo_info "link ${dirinst}/lib/libHepMC3.dylib"
		[ -e ${dirinst}/lib/libHepMC3.dylib ] && ln -s ${dirinst}/lib/libHepMC3.dylib ${dirinst}/lib/libHepMC.dylib
		[ -e ${dirinst}/lib/libHepMC3.so ] && ln -s ${dirinst}/lib/libHepMC3.so ${dirinst}/lib/libHepMC.so
		cd ${cdir}
	fi
fi

if [ -d ${dirinst} ]; then
	export HEPMC3_DIR=${dirinst}
	export HEPMC3_ROOT_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi

cd ${cdir}