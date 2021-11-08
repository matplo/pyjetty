#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os

import pythia8
import pythiaext
import pythiafjext

from heppy.pythiautils import configuration as pyconf

from pyjetty.mputils import mputils

import ROOT
ROOT.gROOT.SetBatch(1)


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

	max_eta_hadron=3
	jet_R0 = 0.4
	# jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(125.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)

	fj.ClusterSequence.print_banner()
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)

	qweights = [1./iw/10. for iw in reversed(range(1, 3)) ]
	qweights.insert(0, 0)

	rout = ROOT.TFile('gen_quench_out.root', 'recreate')
	hpt = []
	for i, w in enumerate(qweights):
		hname = 'hpt_{}'.format(i)
		htitle = 'hpt w={}'.format(w)
		h = ROOT.TH1F(hname, htitle, 10, mputils.logbins(10, 1000, 10))
		hpt.append(h)

	pbar = tqdm.tqdm(range(args.nev))
	for i in pbar:
		if not pythia.next():
			pbar.update(-1)
			continue

		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected)))

		if len(jets_h) < 1:
			continue

		# do your things with jets here...
		for j in jets_h:
			for i, w in enumerate(qweights):
				if w > 0:
					_j = j * w
				else:
					_j = j
				jpt = 0
				print(j.perp(), _j.perp(), w)
				for p in j.constituents():
					if w > 0:
						bp = p.unboost(_j)
					else:
						bp = p
					print(' -', w, p.perp(), bp.perp())
					jpt = jpt + bp.perp()
				hpt[i].Fill(jpt)

	rout.cd()
	rout.Write()
	rout.Close()

	pythia.stat()
	pythia.settings.writeFile(args.py_cmnd_out)


if __name__ == '__main__':
	main()
