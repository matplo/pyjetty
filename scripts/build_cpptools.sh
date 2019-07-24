#!/bin/bash

function os_linux()
{
	_system=$(uname -a | cut -f 1 -d " ")
	if [ $_system == "Linux" ]; then
		echo "yes"
	else
		echo
	fi
}

function os_darwin()
{
	_system=$(uname -a | cut -f 1 -d " ")
	if [ $_system == "Darwin" ]; then
		echo "yes"
	else
		echo
	fi
}

function n_cores()
{
	local _ncores="1"
	[ $(os_darwin) ] && local _ncores=$(system_profiler SPHardwareDataType | grep "Number of Cores" | cut -f 2 -d ":" | sed 's| ||')
	[ $(os_linux) ] && local _ncores=$(lscpu | grep "CPU(s):" | head -n 1 | cut -f 2 -d ":" | sed 's| ||g')
	#[ ${_ncores} -gt "1" ] && retval=$(_ncores-1)
	echo ${_ncores}
}

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

function abspath()
{
  case "${1}" in
    [./]*)
    echo "$(cd ${1%/*}; pwd)/${1##*/}"
    ;;
    *)
    echo "${PWD}/${1}"
    ;;
  esac
}
export -f abspath

function get_opt()
{
    all_opts="$@"
    # echo "options in function: ${all_opts}"
    opt=${1}
    # echo "checking for [${opt}]"
    #opts=("${all_opts[@]:2}")
    opts=$(echo ${all_opts} | cut -d ' ' -f 2-)
    retval=""
    is_set=""
    # echo ".. in [${opts}]"
    for i in ${opts}
    do
    case $i in
        --${opt}=*)
        retval="${i#*=}"
        shift # past argument=value
        ;;
        --${opt})
        is_set=yes
        shift # past argument with no value
        ;;
        *)
            # unknown option
        ;;
    esac
    done
    if [ -z ${retval} ]; then
        echo ${is_set}
    else
        echo ${retval}
    fi
}
export -f get_opt

need_help=$(get_opt "help" $@)
if [ ! -z ${need_help} ]; then
	echo "$0 [--help] [--unsetpyhonpath] [--clean] [--cleanall]"
	exit 0
fi

unsetpyhonpath=$(get_opt "unsetpyhonpath" $@)
if [ ! -z ${unsetpyhonpath} ]; then
    unset PYTHONPATH	&& echo "unsetting PYTHONPATH"
fi

install_path=${SCRIPTPATH}/..
build_path=/tmp/pyjetty_build
echo "[i] building in ${build_path}"

clean=$(get_opt "clean" $@)
if [ ! -z ${clean} ]; then
	echo "[i] removing ${build_path}"
	rm -rf ${build_path}
	exit 0
fi

cleanall=$(get_opt "cleanall" $@)
if [ ! -z ${cleanall} ]; then
	echo "[i] removing ${build_path}"
	rm -rf ${build_path}
	echo "[i] removing ${SCRIPTPATH}/../lib"
	rm -rf ${SCRIPTPATH}/../lib
	exit 0
fi

mkdir -p ${build_path}
if [ -d ${build_path} ]; then
	cd ${build_path}
	cmake -Bbuild -DBUILD_PYTHON=ON -DCMAKE_INSTALL_PREFIX=${install_path} -DCMAKE_BUILD_TYPE=Release $(abspath ${SCRIPTPATH}/../cpptools) \
	&& cmake --build build --target all -- -j $(n_cores) \
	&& cmake --build build --target install
	cd -
else
	echo "[error] unable to access build path: ${build_path}"
fi
