#!/bin/bash

function write_tdraw_file()
{
	inputdir=${1}
	outputfile=${2}
	rm -rf ${outputfile}
	echo "[options]" >> ${outputfile}
	echo "libs =" >> ${outputfile}
	echo "" >> ${outputfile}
	echo "[pt]" >> ${outputfile}
	echo "        active = True" >> ${outputfile}
	echo "        output_file = +_tdraw" >> ${outputfile}
	echo "        input_dir = ${inputdir}/" >> ${outputfile}
	echo "        input_file = pPb_npdf_compare_*.root" >> ${outputfile}
	echo "        tree_name = tnpart" >> ${outputfile}
	echo "        varexp = pt" >> ${outputfile}
	echo "        # selection = eta < -3 && eta > -6" >> ${outputfile}
	#echo "        selection = eta > 3 && eta < 6" >> ${outputfile}
	#echo "        selection = eta > 2.5 && eta < 4." >> ${outputfile}
	echo "        selection = eta < 1." >> ${outputfile}
	echo "        option = e" >> ${outputfile}
	echo "        nentries =" >> ${outputfile}
	echo "        firstentry =" >> ${outputfile}
	echo "        x = 1, 100" >> ${outputfile}
	echo "        nbinsx = 20" >> ${outputfile}
	echo "        logx = True" >> ${outputfile}
	echo "        x_title = "p_{T}"" >> ${outputfile}
	echo "        y_title = counts" >> ${outputfile}
	echo "        title = pt" >> ${outputfile}
	echo "        name = hpt" >> ${outputfile}
}

function write_draw_file()
{
	inputdir=${1}
	outputfile=${2}
	hname=${3}
	rm ${outputfile}
	echo "#-----------------------" >> ${outputfile}
	echo "#figure" >> ${outputfile}
	echo "#geom 500x500" >> ${outputfile}
	echo "#ate" >> ${outputfile}
	echo "#title: smart group" >> ${outputfile}
	echo "#logy 1" >> ${outputfile}
	echo "#logx 1" >> ${outputfile}
	echo "#x p_{T}^{jet} (GeV/c)" >> ${outputfile}
	echo "#y d#sigma/dp_{T} (arb. units)" >> ${outputfile}
	echo "#legend pos=ur title=p-Pb w/ EPPS16 Pb^{208} nPDF" >> ${outputfile}
	echo "#xrange 1, 100" >> ${outputfile}
	echo "##miny 1e-5" >> ${outputfile}
	echo "##maxy 1e3" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_uncerts_output_tdraw.root 		:${hname}_set0_to_graph :p  					: title=nPDF #pi^{0}" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_uncerts_output_tdraw.root		:${hname}_epps16_uncerts :p -k serror noleg +a10 : title=uncerts graph" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_output_tdraw.root 				:${hname} 			:p -k +l9 				: title=free proton PDF #pi^{0}" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_uncerts_output_tdraw.root 	:${hname}_set0_to_graph :p  					: title=nPDF D^{0}" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_uncerts_output_tdraw.root	:${hname}_epps16_uncerts :p -k serror noleg +a10 : title=uncerts graph" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_output_tdraw.root 			:${hname} 			:p -k +l9  				: title=free proton PDF D^{0}" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_uncerts_output_tdraw.root :${hname}_set0_to_graph :p  					: title=nPDF #gamma" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_uncerts_output_tdraw.root	:${hname}_epps16_uncerts :p -k serror noleg +a10 : title=uncerts graph" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_output_tdraw.root 		:${hname} 			:p -k +l9 				: title=free proton PDF #gamma" >> ${outputfile}
	echo
	echo "#-----------------------" >> ${outputfile}
	echo "#figure" >> ${outputfile}
	echo "#geom 500x500" >> ${outputfile}
	echo "#date" >> ${outputfile}
	echo "#title: smart group" >> ${outputfile}
	echo "##logy 1" >> ${outputfile}
	echo "#logx 1" >> ${outputfile}
	echo "#x p_{T}^{jet} (GeV/c)" >> ${outputfile}
	echo "#y ratio  " >> ${outputfile}
	echo "#legend pos=up title=p-Pb EPPS16 Pb^{208} nPDF/ proton PDF (default PYTHIA8.235),, tx_size=0.04" >> ${outputfile}
	echo "##comment item=anti-k_{T} R=0.4 jets |#eta_{jet}| < 3 tx_size=0.04 pos=up" >> ${outputfile}
	echo "#miny 0" >> ${outputfile}
	echo "#maxy 2" >> ${outputfile}
	echo "#xrange 1, 100" >> ${outputfile}
	echo "##line 1, 1, 20, 1, 1, 9, 1" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_uncerts_ratio_tdraw.root         :${hname}_epps16_uncerts_ratio_noErr :p +w2 +k1 : title=#pi^{0} #hat{p_{T}} > 2 " >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_uncerts_ratio_tdraw.root         :${hname}_epps16_uncerts_ratio :p serror noleg -k: title=Graph" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_uncerts_ratio_tdraw.root   :${hname}_epps16_uncerts_ratio_noErr :p +w2 +k2 : title=D^{0} #hat{p_{T}} > 2 " >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_uncerts_ratio_tdraw.root   :${hname}_epps16_uncerts_ratio :p serror noleg -k: title=Graph" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_uncerts_ratio_tdraw.root  :${hname}_epps16_uncerts_ratio_noErr :p +w2 +k4 : title=#gamma #hat{p_{T}} > 2 " >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_uncerts_ratio_tdraw.root  :${hname}_epps16_uncerts_ratio :p serror noleg -k: title=Graph" >> ${outputfile}
	echo "#-----------------------" >> ${outputfile}
	echo "#figure" >> ${outputfile}
	echo "#geom 500x500" >> ${outputfile}
	echo "#date" >> ${outputfile}
	echo "#title: smart group" >> ${outputfile}
	echo "##logy 1" >> ${outputfile}
	echo "#logx 1" >> ${outputfile}
	echo "#x p_{T}^{jet} (GeV/c)" >> ${outputfile}
	echo "#y ratio  " >> ${outputfile}
	echo "#legend pos=up title=p-Pb EPPS16 Pb^{208} nPDF/ proton PDF (default PYTHIA8.235),, tx_size=0.04" >> ${outputfile}
	echo "##comment item=anti-k_{T} R=0.4 jets |#eta_{jet}| < 3 tx_size=0.04 pos=up" >> ${outputfile}
	echo "#miny 0" >> ${outputfile}
	echo "#maxy 2" >> ${outputfile}
	echo "#xrange 1, 100" >> ${outputfile}
	echo "##line 1, 1, 20, 1, 1, 9, 1" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_uncerts_ratio_tdraw.root         :${hname}_epps16_uncerts_ratio_noErr :p +w2 +k1 : title=#pi^{0} #hat{p_{T}} > 2 " >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_uncerts_ratio_tdraw.root         :${hname}_epps16_uncerts_ratio :p serror noleg -k: title=Graph" >> ${outputfile}
	echo "#-----------------------" >> ${outputfile}
	echo "#figure" >> ${outputfile}
	echo "#geom 500x500" >> ${outputfile}
	echo "#date" >> ${outputfile}
	echo "#title: smart group" >> ${outputfile}
	echo "##logy 1" >> ${outputfile}
	echo "#logx 1" >> ${outputfile}
	echo "#x p_{T}^{jet} (GeV/c)" >> ${outputfile}
	echo "#y ratio  " >> ${outputfile}
	echo "#legend pos=up title=p-Pb EPPS16 Pb^{208} nPDF/ proton PDF (default PYTHIA8.235),, tx_size=0.04" >> ${outputfile}
	echo "##comment item=anti-k_{T} R=0.4 jets |#eta_{jet}| < 3 tx_size=0.04 pos=up" >> ${outputfile}
	echo "#miny 0" >> ${outputfile}
	echo "#maxy 2" >> ${outputfile}
	echo "#xrange 1, 100" >> ${outputfile}
	echo "##line 1, 1, 20, 1, 1, 9, 1" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_uncerts_ratio_tdraw.root   :${hname}_epps16_uncerts_ratio_noErr :p +w2 +k2 : title=D^{0} #hat{p_{T}} > 2 " >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_charm_uncerts_ratio_tdraw.root   :${hname}_epps16_uncerts_ratio :p serror noleg -k: title=Graph" >> ${outputfile}
	echo "#-----------------------" >> ${outputfile}
	echo "#figure" >> ${outputfile}
	echo "#geom 500x500" >> ${outputfile}
	echo "#date" >> ${outputfile}
	echo "#title: smart group" >> ${outputfile}
	echo "##logy 1" >> ${outputfile}
	echo "#logx 1" >> ${outputfile}
	echo "#x p_{T}^{jet} (GeV/c)" >> ${outputfile}
	echo "#y ratio  " >> ${outputfile}
	echo "#legend pos=up title=p-Pb EPPS16 Pb^{208} nPDF/ proton PDF (default PYTHIA8.235),, tx_size=0.04" >> ${outputfile}
	echo "##comment item=anti-k_{T} R=0.4 jets |#eta_{jet}| < 3 tx_size=0.04 pos=up" >> ${outputfile}
	echo "#miny 0" >> ${outputfile}
	echo "#maxy 2" >> ${outputfile}
	echo "#xrange 1, 100" >> ${outputfile}
	echo "##line 1, 1, 20, 1, 1, 9, 1" >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_uncerts_ratio_tdraw.root  :${hname}_epps16_uncerts_ratio_noErr :p +w2 +k4 : title=#gamma #hat{p_{T}} > 2 " >> ${outputfile}
	echo "${inputdir}/pPb_npdf_compare_photon_uncerts_ratio_tdraw.root  :${hname}_epps16_uncerts_ratio :p serror noleg -k: title=Graph" >> ${outputfile}
}

indir=${1}

[ -z ${indir} ] && echo "usage: $(basename ${0}) <input_directory>" && exit 0
[ ! -d ${indir} ] && echo "usage: $(basename ${0}) <input_directory>" && exit 0

tdraw_file=$(basename ${indir})
tdraw_file="tdraw_${tdraw_file}.cfg"
write_tdraw_file ${indir} ${tdraw_file}
tdraw_cfg.py ${tdraw_file} --clean

for opt in "_charm" "_photon" ""
do
	./make_uncert_epps16_tdraw.py -i ${indir}/pPb_npdf_compare${opt}_EPPS16_set0_output_tdraw.root
done

draw_file=$(basename ${indir})
draw_file="${tdraw_file}.draw"
write_draw_file ${indir} ${draw_file} pt

echo "[i] ${draw_file} written."
