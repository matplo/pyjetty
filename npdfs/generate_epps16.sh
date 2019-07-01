#!/bin/bash

outdir=${PWD}/epps16-pp-5000-bias/
mkdir -p ${outdir}
./jets_epps16.py -w ${outdir}/jets_npdf_compare.dat --allset --ecm 5000. --biaspow 4 --biasref 20
