#!/bin/bash

srcdir=${1}
wdir=${2}
if [ ! -z ${wdir} ]; then
	[ ! -d ${wdir} ] && mkdir -pv ${wdir}
fi

if [ -d ${srcdir} ]; then
	if [ -d ${wdir} ]; then
		cd ${wdir}
		[ ! -e fjcontrib-1.041.tar.gz ] && wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.041.tar.gz
		if [ -e fjcontrib-1.041.tar.gz ]; then
			if [ ! -d RecursiveTools ]; then
				tar zxvf fjcontrib-1.041.tar.gz RecursiveTools
				patch ${srcdir}/RecursiveTools/RecursiveSymmetryCutBase.hh -i ${srcdir}/patches/RecursiveSymmetryCutBase.patch
			fi

			if [ ! -d LundPlane ]; then
				tar zxvf fjcontrib-1.041.tar.gz LundPlane
				rm ${srcdir}/LundPlane/example_*.cc
				patch ${srcdir}/LundPlane/SecondaryLund.hh -i ${srcdir}/patches/SecondaryLund.patch
				patch ${srcdir}/LundPlane/LundGenerator.hh -i ${srcdir}/patches/LundGenerator.patch
			fi
		fi
	fi
fi
