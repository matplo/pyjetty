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
[ -z ${version} ] && version=8235
warning "pythia8 version ${version}"
fname=pythia${version}
dirsrc=${SCRIPTPATH}/build/pythia${version}
dirinst=${SCRIPTPATH}/packages/pythia${version}

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

if [ ! -e ${SCRIPTPATH}/build/${fname}.tgz ]; then
	cd ${SCRIPTPATH}/build
	wget http://home.thep.lu.se/~torbjorn/pythia8/${fname}.tgz
fi

if [ ! -d ${dirsrc} ]; then
	cd ${SCRIPTPATH}/build
	tar zxvf ${fname}.tgz
fi

function run_configure()
{
	[ ! -z ${LHAPDF6_DIR} ] && lhapd6_opt="	--with-lhapdf6=${LHAPDF6_DIR} --with-lhapdf6-include=${LHAPDF6_DIR}/include --with-lhapdf6-lib=${LHAPDF6_DIR}/lib"
	[ ! -z ${HEPMC2_DIR} ] && hepmc2_opt=" --with-hepmc2=${HEPMC2_DIR} --with-hepmc2-include=${HEPMC2_DIR}/include --with-hepmc2-lib=${HEPMC2_DIR}/lib"
	[ ! -z ${HEPMC3_DIR} ] && hepmc3_opt=" --with-hepmc3=${HEPMC3_DIR} --with-hepmc3-include=${HEPMC3_DIR}/include --with-hepmc3-lib=${HEPMC3_DIR}/lib"
		./configure --prefix=${dirinst} \
			--with-python-include=${python_inc_dir} \
			--with-python-bin=${python_bin_dir} \
			${lhapd6_opt} ${hepmc2_opt} ${hepmc3_opt}
}

redo=$(get_opt "re" $@)
if [ ! -d ${dirinst} ] || [ "x${redo}" == "xyes" ]; then
	if [ -d ${dirsrc} ]; then
		cd ${dirsrc}
	    # echo "unsetting PYTHONPATH"
		[ "$(get_opt "unset" $@)" == "xyes" ] && unset PYTHONPATH && warning "unsetting PYTHONPATH"
	    python_inc_dir=$(python3-config --includes | cut -d' ' -f 1 | cut -dI -f 2)
	    python_exec=$(which python3)
	    python_bin_dir=$(dirname ${python_exec})
	    echo_info "python exec: ${python_exec}"
	    echo_info "python bin dir: ${python_bin_dir}"
	    echo_info "python include: ${python_inc_dir}"
	    if [ ! -e ${python_bin_dir}/python ]; then
	    	python_bin_dir=${SCRIPTPATH}/build/pythia-python-bin
	    	mkdir ${python_bin_dir}
	    	ln -s ${python_exec} ${python_bin_dir}/python 
		    warning "fix-up-python bin dir: ${python_bin_dir}"
		fi
    	if [ ! -d ${dirinst} ]; then
			run_configure
			configure_only=$(get_opt "configure-only" $@)
			[ "x${configure_only}" == "xyes" ] && grace_return && return 0
			make -j $(n_cores) && make install
		fi
		cd - 2>&1 > /dev/null
	fi
fi

if [ -d ${dirinst} ]; then
	export PYTHIA_DIR=${dirinst}
	export PYTHIA8_ROOT_DIR=${dirinst}
	export PATH=$PATH:${dirinst}/bin
	export PYTHONPATH=${PYTHONPATH}:${dirinst}/lib
fi

cd ${cdir}