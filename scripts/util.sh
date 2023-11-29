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

function abspath_python_expand()
{
	rv=$(python -c "import os; print(os.path.abspath(os.path.expandvars(\"${1}\")))")
	echo ${rv}
}
export -f abspath_python_expand

function os_linux()
{
	_system=$(uname -a | cut -f 1 -d " ")
	if [ $_system == "Linux" ]; then
		echo "yes"
	else
		echo
	fi
}
export -f os_linux

function os_darwin()
{
	_system=$(uname -a | cut -f 1 -d " ")
	if [ $_system == "Darwin" ]; then
		echo "yes"
	else
		echo
	fi
}
export -f os_darwin

function n_cores()
{
	local _ncores="1"
	[ $(os_darwin) ] && local _ncores=$(system_profiler SPHardwareDataType | grep "Number of Cores" | cut -f 2 -d ":" | sed 's| ||' | cut -f 1 -d " ")
	[ $(os_linux) ] && local _ncores=$(lscpu | grep "CPU(s):" | head -n 1 | cut -f 2 -d ":" | sed 's| ||g')
	#[ ${_ncores} -gt "1" ] && retval=$(_ncores-1)
	echo ${_ncores}
}
export -f n_cores

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
if [ "x${need_help}" == "xyes" ]; then
    echo "[i] help requested"
fi

# asetting=$(get_opt "asetting" $@)
# if [ ! -z ${asetting} ]; then
#     echo "[i] asetting: ${asetting}"
# fi

function echo_info()
{
	(>&2 echo "[info] $@")
}
export -f echo_info

function echo_warning()
{
	(>&2 echo -e "\033[1;93m$@ \033[0m")
}
export -f echo_warning

function echo_error()
{
	(>&2 echo -e "\033[1;31m$@ \033[0m")
}
export echo_error

function echo_note_red()
{
	(>&2 echo -e "\033[1;31m[note] $@ \033[0m")
}
export echo_note_red

function note_red()
{
	(>&2 echo -e "\033[1;31m[note] $@ \033[0m")
}
export -f note_red

function separator()
{
	echo -e "\033[1;32m$(padding "[ ${1} ]" "-" 50 center) \033[0m"
	## colors at http://misc.flogisoft.com/bash/tip_colors_and_formatting
}
export -f separator

function echo_note()
{
	echo_warning "$(padding "[note] ${@}" "-" 10 left)"
}
export -f echo_note

function note()
{
	echo_warning "$(padding "[note] ${@}" "-" 10 left)"
}
export -f note

function warning()
{
	echo_warning "[warning] $(padding "[${@}] " "-" 40 right)"
}
export -f warning

function error()
{
	echo_error "[error] $(padding "[${@}] " "-" 42 right)"
}
export -f error

function padding ()
{
	CONTENT="${1}";
	PADDING="${2}";
	LENGTH="${3}";
	TRG_EDGE="${4}";
	case "${TRG_EDGE}" in
		left) echo ${CONTENT} | sed -e :a -e 's/^.\{1,'${LENGTH}'\}$/&\'${PADDING}'/;ta'; ;;
		right) echo ${CONTENT} | sed -e :a -e 's/^.\{1,'${LENGTH}'\}$/\'${PADDING}'&/;ta'; ;;
		center) echo ${CONTENT} | sed -e :a -e 's/^.\{1,'${LENGTH}'\}$/'${PADDING}'&'${PADDING}'/;ta'
	esac
	return ${RET__DONE};
}
export -f padding

function setup_python_env()
{
	separator "setup_python_env()"
	[ -z ${HEPPY_USER_PYTHON_VERSION} ] && export HEPPY_USER_PYTHON_VERSION=python
	_PYTHON_EXECUTABLE=$(which ${HEPPY_USER_PYTHON_VERSION})
	if [ -z ${_PYTHON_EXECUTABLE} ]; then
		_PYTHON_VERSION=""
		_PYTHON_BIN_DIR=""
	    _PYTHON_INCLUDE_DIR=""
	    _PYTHON_LIBDIR=""
    	_PYTHON_NUMPY_INCLUDE_DIR=""
	    _PYTHON_CONFIG_EXECUTABLE=""
	    _PYTHON_LIBS=""
	    _PYTHON_CONFIG_LDFLAGS=""
	    _PYTHON_CONFIG_INCLUDES=""
    	_PYTHON_SETUP=""
    	_PYTHON_LIBS_LINK=""
	else
		_PYTHON_VERSION=$(${_PYTHON_EXECUTABLE} --version 2>&1 | cut -f 2 -d' ' | cut -f 1-2 -d.)
		_PYTHON_BIN_DIR=$(dirname ${_PYTHON_EXECUTABLE})
	    _PYTHON_INCLUDE_DIR=$(${_PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")
	    _PYTHON_LIBDIR=$(${_PYTHON_EXECUTABLE} -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
    	_PYTHON_NUMPY_INCLUDE_DIR=$(${_PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())")

	    _PYTHON_CONFIG_EXECUTABLE=$(which ${HEPPY_USER_PYTHON_VERSION}-config)
		if [ -z ${_PYTHON_CONFIG_EXECUTABLE} ]; then
			warning "guessing python libs..."
		    _PYTHON_LIBS="-lpython${_PYTHON_VERSION}"
		else
		    _PYTHON_LIBS=$(${_PYTHON_CONFIG_EXECUTABLE} --libs)
	    	_PYTHON_LIBS_LINK="-L${_PYTHON_LIBDIR} ${_PYTHON_LIBS}"
	    	_PYTHON_CONFIG_LDFLAGS=$(${_PYTHON_CONFIG_EXECUTABLE} --ldflags)
	    	_PYTHON_CONFIG_INCLUDES=$(${_PYTHON_CONFIG_EXECUTABLE} --includes)
	    	_PYTHON_SETUP=TRUE
	    fi
	fi

	export HEPPY_PYTHON_EXECUTABLE=${_PYTHON_EXECUTABLE}
	export HEPPY_PYTHON_VERSION=${_PYTHON_VERSION}
	export HEPPY_PYTHON_BIN_DIR=${_PYTHON_BIN_DIR}
	export HEPPY_PYTHON_CONFIG_INCLUDES=${_PYTHON_CONFIG_INCLUDES}
	export HEPPY_PYTHON_CONFIG_LDFLAGS=${_PYTHON_CONFIG_LDFLAGS}
	export HEPPY_PYTHON_LIBS=${_PYTHON_LIBS}
	export HEPPY_PYTHON_LIBS_LINK=${_PYTHON_LIBS_LINK}
	export HEPPY_PYTHON_CONFIG_EXECUTABLE=${_PYTHON_CONFIG_EXECUTABLE}
	export HEPPY_PYTHON_NUMPY_INCLUDE_DIR=${_PYTHON_NUMPY_INCLUDE_DIR}
	export HEPPY_PYTHON_LIBDIR=${_PYTHON_LIBDIR}
	export HEPPY_PYTHON_INCLUDE_DIR=${_PYTHON_INCLUDE_DIR}
	export HEPPY_PYTHON_SETUP=${_PYTHON_SETUP}
}
export -f setup_python_env

function echo_python_setup()
{
	echo_info "HEPPY_PYTHON_SETUP=${HEPPY_PYTHON_SETUP}"
	echo_info "HEPPY_PYTHON_EXECUTABLE=${HEPPY_PYTHON_EXECUTABLE}"
	echo_info "HEPPY_PYTHON_VERSION=${HEPPY_PYTHON_VERSION}"
	echo_info "HEPPY_PYTHON_BIN_DIR=${HEPPY_PYTHON_BIN_DIR}"
	echo_info "HEPPY_PYTHON_INCLUDE_DIR=${HEPPY_PYTHON_INCLUDE_DIR}"
	echo_info "HEPPY_PYTHON_LIBDIR=${HEPPY_PYTHON_LIBDIR}"
	echo_info "HEPPY_PYTHON_NUMPY_INCLUDE_DIR=${HEPPY_PYTHON_NUMPY_INCLUDE_DIR}"
	echo_info "HEPPY_PYTHON_CONFIG_EXECUTABLE=${HEPPY_PYTHON_CONFIG_EXECUTABLE}"
	echo_info "HEPPY_PYTHON_LIBS=${HEPPY_PYTHON_LIBS}"
	echo_info "HEPPY_PYTHON_CONFIG_LDFLAGS=${HEPPY_PYTHON_CONFIG_LDFLAGS}"
	echo_info "HEPPY_PYTHON_CONFIG_INCLUDES=${HEPPY_PYTHON_CONFIG_INCLUDES}"
	echo_info "HEPPY_PYTHON_LIBS_LINK=${HEPPY_PYTHON_LIBS_LINK}"
}
export -f echo_python_setup

function add_path()
{
	path=${1}
	if [ ! -z ${path} ] && [ -d ${path} ]; then
		echo_info "adding ${path} to PATH"
		if [ -z ${PATH} ]; then
			export PATH=${path}
		else
			export PATH=${path}:${PATH}
		fi
	else
		[ "x${2}" == "xdebug" ] && echo_error "ignoring ${path} for PATH"
	fi
}
export -f add_path

function add_pythonpath()
{
	path=${1}
	if [ ! -z ${path} ] && [ -d ${path} ]; then
		echo_info "adding ${path} PYTHONPATH"
		if [ -z ${PYTHONPATH} ]; then
			export PYTHONPATH=${path}
		else
			export PYTHONPATH=${path}:${PYTHONPATH}
		fi
	else
		[ "x${2}" == "xdebug" ] && echo_error "ignoring ${path} for PYTHONPATH"
	fi
}
export -f add_pythonpath

function add_ldpath()
{
	path=${1}
	if [ ! -z ${path} ] && [ -d ${path} ]; then
		echo_info "adding ${path} to LD_LIBRARY_PATH"
		if [ -z ${LD_LIBRARY_PATH} ]; then
			export LD_LIBRARY_PATH=${path}
		else
			export LD_LIBRARY_PATH=${path}:${LD_LIBRARY_PATH}
		fi
	else
		[ "x${2}" == "xdebug" ] && echo_error "ignoring ${path} for LD_LIBRARY_PATH"
	fi
}
export -f add_ldpath

function add_dyldpath()
{
	path=${1}
	if [ ! -z ${path} ] && [ -d ${path} ]; then
		echo_info "adding ${path} to DYLD_LIBRARY_PATH"
		if [ -z ${DYLD_LIBRARY_PATH} ]; then
			export DYLD_LIBRARY_PATH=${path}
		else
			export DYLD_LIBRARY_PATH=${path}:${DYLD_LIBRARY_PATH}
		fi
	else
		[ "x${2}" == "xdebug" ] && echo_error "ignoring ${path} for DYLD_LIBRARY_PATH"
	fi
}
export -f add_dyldpath

function add_path_module()
{
	path=${@:3}
	what=${2}
	modulefile=${1}
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
		[ ${HEPPY_DEBUG} ] && warning "add_path_module:: ignoring ${path} for ${what}"
	fi
}
export -f add_path_module

function setenv_module()
{
	path=${@:3}
	what=${2}
	modulefile=${1}
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
			echo "setenv ${what} \"${path}\"" >> ${modulefile}
		fi
	else
		[ ${HEPPY_DEBUG} ] && warning "setenv_module:: ignoring ${path} for ${what}"
	fi
}
export -f setenv_module

function setalias_module()
{
	path=${@:3}
	what=${2}
	modulefile=${1}
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
