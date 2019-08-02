if [ ! -z ${PYJETTY_EXTERNAL_UTIL_LOADED} ]; then
	return 0
fi
PYJETTY_EXTERNAL_UTIL_LOADED=yes

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
	_PYTHON_EXECUTABLE=$(which python)
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
		_PYTHON_VERSION=$(${_PYTHON_EXECUTABLE} --version | cut -f 2 -d' ' | cut -f 1-2 -d.)
		_PYTHON_BIN_DIR=$(dirname ${_PYTHON_EXECUTABLE})
	    _PYTHON_INCLUDE_DIR=$(${_PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")
	    _PYTHON_LIBDIR=$(${_PYTHON_EXECUTABLE} -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
    	_PYTHON_NUMPY_INCLUDE_DIR=$(${_PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())")

	    _PYTHON_CONFIG_EXECUTABLE=$(which python-config)
		if [ -z ${_PYTHON_CONFIG_EXECUTABLE} ]; then
			warning "guessing python libs..."
		    _PYTHON_LIBS="-lpython${_PYTHON_VERSION}"
		else
		    _PYTHON_LIBS=$(${_PYTHON_CONFIG_EXECUTABLE} --libs)
	    	_PYTHON_LIBS_LINK="-L${_PYTHON_LIBDIR} ${_PYTHON_LIBS}"
	    	_PYTHON_CONFIG_LDFLAGS=$(${_PYTHON_CONFIG_EXECUTABLE} --ldflags)
	    	_PYTHON_CONFIG_INCLUDES=$(python-config --includes)
	    	_PYTHON_SETUP=TRUE
	    fi
	fi

	export PYJETTY_PYTHON_EXECUTABLE=${_PYTHON_EXECUTABLE}
	export PYJETTY_PYTHON_VERSION=${_PYTHON_VERSION}
	export PYJETTY_PYTHON_BIN_DIR=${_PYTHON_BIN_DIR}
	export PYJETTY_PYTHON_CONFIG_INCLUDES=${_PYTHON_CONFIG_INCLUDES}
	export PYJETTY_PYTHON_CONFIG_LDFLAGS=${_PYTHON_CONFIG_LDFLAGS}
	export PYJETTY_PYTHON_LIBS=${_PYTHON_LIBS}
	export PYJETTY_PYTHON_LIBS_LINK=${_PYTHON_LIBS_LINK}
	export PYJETTY_PYTHON_CONFIG_EXECUTABLE=${_PYTHON_CONFIG_EXECUTABLE}
	export PYJETTY_PYTHON_NUMPY_INCLUDE_DIR=${_PYTHON_NUMPY_INCLUDE_DIR}
	export PYJETTY_PYTHON_LIBDIR=${_PYTHON_LIBDIR}
	export PYJETTY_PYTHON_INCLUDE_DIR=${_PYTHON_INCLUDE_DIR}
	export PYJETTY_PYTHON_SETUP=${_PYTHON_SETUP}
}
export -f setup_python_env

function echo_python_setup()
{
	echo_info "PYJETTY_PYTHON_SETUP=${PYJETTY_PYTHON_SETUP}"
	echo_info "PYJETTY_PYTHON_EXECUTABLE=${PYJETTY_PYTHON_EXECUTABLE}"
	echo_info "PYJETTY_PYTHON_VERSION=${PYJETTY_PYTHON_VERSION}"
	echo_info "PYJETTY_PYTHON_BIN_DIR=${PYJETTY_PYTHON_BIN_DIR}"
	echo_info "PYJETTY_PYTHON_INCLUDE_DIR=${PYJETTY_PYTHON_INCLUDE_DIR}"
	echo_info "PYJETTY_PYTHON_LIBDIR=${PYJETTY_PYTHON_LIBDIR}"
	echo_info "PYJETTY_PYTHON_NUMPY_INCLUDE_DIR=${PYJETTY_PYTHON_NUMPY_INCLUDE_DIR}"
	echo_info "PYJETTY_PYTHON_CONFIG_EXECUTABLE=${PYJETTY_PYTHON_CONFIG_EXECUTABLE}"
	echo_info "PYJETTY_PYTHON_LIBS=${PYJETTY_PYTHON_LIBS}"
	echo_info "PYJETTY_PYTHON_CONFIG_LDFLAGS=${PYJETTY_PYTHON_CONFIG_LDFLAGS}"
	echo_info "PYJETTY_PYTHON_CONFIG_INCLUDES=${PYJETTY_PYTHON_CONFIG_INCLUDES}"
	echo_info "PYJETTY_PYTHON_LIBS_LINK=${PYJETTY_PYTHON_LIBS_LINK}"
}
export -f echo_python_setup
