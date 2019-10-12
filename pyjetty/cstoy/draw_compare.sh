#!/bin/bash

dfile="tmp.draw"
echo "#figure" | tee ${dfile}
echo "#legend pos=ur alpha=50 title=<lt> (lines: matched)" | tee -a ${dfile}
echo "#xrange 0, 1.2" | tee -a ${dfile}
echo "<f1> : <h1> : p EX0 +m75: title=<t1>" | tee -a ${dfile}
echo "<f1> : <h2> : p EX0 +m75: title=<t2>" | tee -a ${dfile}
echo "<f1> : <h2x> : hist l -k -l noleg : title=<t2>" | tee -a ${dfile}
echo "<f2> : <h2> : p EX0 +m75: title=<t3>" | tee -a ${dfile}
echo "<f2> : <h2x> : hist l -k -l noleg : title=<t3>" | tee -a ${dfile}
echo "<f3> : <h2> : p EX0 +m75: title=<t4>" | tee -a ${dfile}
echo "<f3> : <h2x> : hist l -k -l noleg : title=<t4>" | tee -a ${dfile}

for pt in 20 40 60 80 100
do
	if [ "${pt}" == "20" ]; then
		extra="-r pos=ur:pos=ul"
	else
		extra=""
	fi
	[ "${1}" == "gui" ] && gui_draw_select.py ${dfile} --wname draw_compare_${pt} \
		-r lt:"${pt}<p_{T}<${pt}+20" \
		-r f1:output_alpha_0_dRmax_0.2_SDzcut_0.2_emb_merged_hout_thg.root \
		-r h1:h_thetagv_pt${pt} \
		-r h2:h_thetag_pt${pt}  \
		-r h2x:h_thetag_pt${pt}_match  \
		-r t1:vacuum \
		-r t2:dR_{max}=0.2 \
		-r t3:dR_{max}=0.25 \
		-r t4:dR_{max}=0.4 \
		-r f2:output_alpha_0_dRmax_0.25_SDzcut_0.2_emb_merged_hout_thg.root \
		-r f3:output_alpha_0_dRmax_0.4_SDzcut_0.2_emb_merged_hout_thg.root 	\
		${extra} \
		--preent pdf1 --quit
done
