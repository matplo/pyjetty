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
THISD=$(thisdir)
source ${THISD}/../scripts/util.sh
separator "heppy: ${BASH_SOURCE}"

install_path=$(abspath ${THISD})
build_path=${THISD}/build

clean=$(get_opt "clean" $@)
if [ ! -z ${clean} ]; then
    rm -rf ${build_path} ${THISD}/lib ${THISD}/lib64 
fi

verbose=$(get_opt "verbose" $@)
if [ ! -z ${verbose} ]; then
    export VERBOSE=1
else
    unset VERBOSE
fi

build_python_iface="TRUE"
cxx_only=$(get_opt "cxx-only" $@)
if [ "x${cxx_only}" == "xyes" ]; then
    build_python_iface="FALSE"
fi

build_configuration="Release"
debug=$(get_opt "debug" $@)
if [ ! -z ${debug} ]; then
    build_configuration="Debug"
fi

build_tenngen=""
tengen=$(get_opt "tenngen" $@)
if [ ! -z ${tengen} ]; then
    build_tenngen="-DTENNGEN=TRUE"
    echo_note "TennGen will be build"
    if [ ! -d "${THISD}/src/TennGen" ]; then
        cdir=$PWD
        cd ${THISD}/src/
        git clone https://github.com/matplo/TennGen.git
        cd ${cdir}
    fi
else
    echo_note "TennGen will NOT be build"    
fi

build_tglaubermc=""
tengen=$(get_opt "tglaubermc" $@)
if [ ! -z ${tengen} ]; then
    build_tglaubermc="-DTGLAUBERMC=TRUE"
    echo_note "TGlauberMC will be build"
    if [ ! -d "${THISD}/src/TGlauberMC" ]; then
        cdir=$PWD
        cd ${THISD}/src/
        git clone https://github.com/matplo/TGlauberMC.git
        cd ${cdir}
    fi
else
    echo_note "TGlauberMC will NOT be build"    
fi

configure_only=$(get_opt "configure-only" $@)

RFEAT=$(root-config --features)
echo $RFEAT
cxxstd=cxx14
if grep -q "$cxxstd" <<< "$RFEAT"; then
    CXX_STD=-DCMAKE_CXX_STANDARD=14
    warning "setting $CXX_STD because or ROOT w/ cxx14"
fi

echo "[i] building in ${build_path}"
mkdir -p ${build_path}
if [ -d ${build_path} ]; then
    cd ${build_path}
    separator "configure"
    cmake -DMAKE_MODULE=TRUE -DBUILD_PYTHON=${build_python_iface} ${build_python_iface} \
            -DCMAKE_INSTALL_PREFIX=${install_path} -DCMAKE_BUILD_TYPE=${build_configuration} \
            ${build_tenngen} \
            ${build_tglaubermc} \
	    ${CXX_STD} \
            ${THISD}
    if [ "x${configure_only}" == "xyes" ]; then
        warning "stopping short of building...- configure-only requested"
        exit 0
    fi
    separator "build"
    cmake  --build . --target all -- -j $(n_cores) && cmake --build . --target install
    # separator "install"
else
	error "unable to access build path: ${build_path}"
fi
