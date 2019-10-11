#!/bin/bash

files=$(ls *_hout_*.root)
for fn in ${files}
do
	echo ${fn}
	gui_draw_select.py cstoy_std_file.draw -r f1:${fn} --wname ${fn} &
done
