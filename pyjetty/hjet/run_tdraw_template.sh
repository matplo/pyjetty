#!/bin/bash

pthats="4 11 21 36 56 84 117 156 200 249 1000"

for pthat in ${pthats}
do
	merged=hout_pthat_${pthat}.root
	[ "x${1}" == "xoverwrite" ] && rm -vf ${merged}
	[ ! -f ${merged} ] && tdraw_cfg.py tdraw_hjet_template.cfg -r pthat:${pthat} --clean --merged-output ${merged}
done

rm -vf hout_pthat_merged.root
hadd -f hout_pthat_merged.root hout_pthat_*.root
