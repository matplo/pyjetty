#!/bin/bash

dfile="tmp.draw"
echo "#figure" | tee ${dfile}
echo "#legend pos=ur alpha=0 title=<lt>" | tee -a ${dfile}
echo "#xrange 0, 1.2" | tee -a ${dfile}
echo "<f1> : <h1> : p EX0 : title=<t1>" | tee -a ${dfile}
echo "<f2> : <h2> : p EX0 : title=<t2>" | tee -a ${dfile}
echo "<f3> : <h3> : p EX0 : title=<t3>" | tee -a ${dfile}
echo "<f4> : <h4> : p EX0 : title=<t4>" | tee -a ${dfile}

for pt in 20 40 60 80
do
	if [ "${pt}" == "20" ]; then
		extra="-r pos=ur:pos=lr"
	else
		extra=""
	fi
	[ "${1}" == "gui" ] && gui_draw_select.py ${dfile} --wname draw_compare_${pt} \
		-r lt:"${pt}<p_{T}<${pt}+20" \
		-r f1:output_alpha_0_dRmax_0.2_SDzcut_0.2_emb_merged_hout_thg.root 	-r h1:h_thetagv_pt${pt} 	-r t1:vacuum \
		-r f2:output_alpha_0_dRmax_0.2_SDzcut_0.2_emb_merged_hout_thg.root 	-r h2:h_thetag_pt${pt} 		-r t2:dR_{max}=0.2 \
		-r f3:output_alpha_0_dRmax_0.25_SDzcut_0.2_emb_merged_hout_thg.root -r h3:h_thetag_pt${pt} 		-r t3:dR_{max}=0.25 \
		-r f4:output_alpha_0_dRmax_0.4_SDzcut_0.2_emb_merged_hout_thg.root 	-r h4:h_thetag_pt${pt} 		-r t4:dR_{max}=0.4 \
		${extra} \
		--preent pdf1 --quit
done
