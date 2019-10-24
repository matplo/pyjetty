#!/bin/bash

indir=/Users/ploskon/devel/pyjetty/pyjetty/hjet/hiccupgen/hjet/pthat
files=$(find ${indir} -name "h_jet_ch_R04_tranges_6-7_20-30_runid_*_pthatmin_*_hard.root")
echo ${files} | wc -w
parallel ./analyze_hjet.py -t {1} -f {2} ::: hjetT_6_7 hjetT_20_30 evT ::: ${files}
