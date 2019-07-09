#!/bin/bash

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

function usage()
{
	echo "usage:" $(basename ${0}) "--mode=bias|pthat [--power=<power>] [--ref=<ref>] [--pthatmin=<minval>] [--nev=<n_events>]"
}

smode=$(get_opt "mode" $@)
biaspower=$(get_opt "power" $@)
biasref=$(get_opt "ref" $@)
pthatmin=$(get_opt "pthatmin" $@)
nev=$(get_opt "nev" $@)
if [ -z ${nev} ]; then
    nev=1000
    echo "[i] setting default nev: ${nev}"
fi

if [ -z ${smode} ]; then
	usage
	exit 0
else
	[ "x${smode}" == "xbias" ]  && ( [ -z ${biaspower} ] || [ -z ${biasref} ] ) && usage	&& exit 0
	[ "x${smode}" == "xpthat" ] && [ -z ${pthatmin} ] && usage && exit 0
fi


if [ "x${smode}" == "xbias" ]; then
	[ -z ${biaspower} ] && biaspower=4
	[ -z ${biasref} ] && biasref=5
	echo "mode : ${smode} && biaspower : ${biaspower} && biasref : ${biasref}"
	outdir="${PWD}/epps16-PbPb-ehigh-bias-${biaspower}-${biasref}"
	mkdir -p ${outdir}
	sets=$(seq -1 0)
	sets=$(seq -1 40)
	parallel --bar ./PbPb_epps16.py -g ${outdir}/PbPb_npdf_compare.dat --epps16set {} --ecm high --biaspow 4 --biasref 20 -n ${nev} ::: ${sets}
	parallel --bar ./PbPb_epps16.py -g ${outdir}/PbPb_npdf_compare_charm.dat --charm --epps16set {} --ecm high --biaspow 4 --biasref 20 -n ${nev} ::: ${sets}
	parallel --bar ./PbPb_epps16.py -g ${outdir}/PbPb_npdf_compare_photon.dat --photon --epps16set {} --ecm high --biaspow 4 --biasref 20 -n ${nev} ::: ${sets}
else
	[ -z ${pthatmin} ] && pthatmin=10
	echo "mode : ${mode} && pthatmin : ${pthatmin}"
	if [ ! -z ${pthatmin} ]; then
		outdir=${PWD}/epps16-PbPb-ehigh-pthat${pthatmin}/
		mkdir -p ${outdir}
		sets=$(seq -1 0)
		sets=$(seq -1 40)
		parallel --bar ./PbPb_epps16.py -g ${outdir}/PbPb_npdf_compare.dat --epps16set {} --ecm high --pthatmin ${pthatmin} -n ${nev} ::: ${sets}
		parallel --bar ./PbPb_epps16.py -g ${outdir}/PbPb_npdf_compare_charm.dat --charm --epps16set {} --ecm high --pthatmin ${pthatmin} -n ${nev} ::: ${sets}
		parallel --bar ./PbPb_epps16.py -g ${outdir}/PbPb_npdf_compare_photon.dat --photon --epps16set {} --ecm high --pthatmin ${pthatmin} -n ${nev} ::: ${sets}
	fi
fi
