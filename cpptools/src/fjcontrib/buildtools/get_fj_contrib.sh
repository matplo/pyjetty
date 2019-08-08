#!/bin/bash

srcdir=${1}
wdir=${2}
if [ ! -z ${wdir} ]; then
	[ ! -d ${wdir} ] && mkdir -p ${wdir}
fi

if [ -d ${srcdir} ]; then
	if [ -d ${wdir} ]; then
		cd ${wdir}
		[ ! -e fjcontrib-1.041.tar.gz ] && wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.041.tar.gz
		cd -
		if [ -e ${wdir}/fjcontrib-1.041.tar.gz ]; then
			if [ ! -d ${srcdir}/fjcontrib-1.041/RecursiveTools ]; then
				cd ${srcdir}
				tar zxvf ${wdir}/fjcontrib-1.041.tar.gz fjcontrib-1.041/RecursiveTools
				rm fjcontrib-1.041/RecursiveTools/example_*.cc
				patch fjcontrib-1.041/RecursiveTools/RecursiveSymmetryCutBase.hh -i ${srcdir}/patches/RecursiveSymmetryCutBase.patch
				cp -v ${srcdir}/custom/* fjcontrib-1.041/RecursiveTools/
			fi

			if [ ! -d ${srcdir}/fjcontrib-1.041/LundPlane ]; then
				cd ${srcdir}
				tar zxvf ${wdir}/fjcontrib-1.041.tar.gz fjcontrib-1.041/LundPlane
				rm fjcontrib-1.041/LundPlane/example_*.cc
				patch fjcontrib-1.041/LundPlane/SecondaryLund.hh -i ${srcdir}/patches/SecondaryLund.patch
				patch fjcontrib-1.041/LundPlane/LundGenerator.hh -i ${srcdir}/patches/LundGenerator.patch
			fi
		fi
	fi
fi
