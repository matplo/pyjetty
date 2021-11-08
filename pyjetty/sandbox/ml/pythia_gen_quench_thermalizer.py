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
ROOT.gSystem.Load('libpyjetty_rutil')


def main():
	parser = argparse.ArgumentParser(description='pythia8 in python', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	

	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

	max_eta_hadron=3
	jet_R0 = 0.4
	jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(125.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	# jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(max_eta_hadron - 1.05 * jet_R0)
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)

	fj.ClusterSequence.print_banner()
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)

	qweights = [2, 4, 6, 10] # n-parts quenched
	qweights.insert(0, 0)

	rout = ROOT.TFile('gen_quench_out.root', 'recreate')
	hpt = []
	hz = []
	thm = [] 
	tn = []
	for i, w in enumerate(qweights):
		hname = 'hpt_{}'.format(i)
		htitle = 'hpt w={}'.format(w)
		rout.cd()
		# h = ROOT.TH1F(hname, htitle, 10, mputils.logbins(10, 500, 10))
		h = ROOT.TH1F(hname, htitle, 10, 0, 250)
		hpt.append(h)
		hname = 'hpz_{}'.format(i)
		htitle = 'hpz w={}'.format(w)
		h = ROOT.TH1F(hname, htitle, 10, 0, 1)
		hz.append(h)
		hname = 'th_{}'.format(i)
		htitle = 'thf w={}'.format(w)
		rout.cd()
		_t = ROOT.RUtil.Thermalizer(0.7, w, 1.0, max_eta_hadron)
		thm.append(_t)

		_tn = ROOT.TNtuple('tree_Particle_gen_{}'.format(w), 'particles from thermalizer {}'.format(w), 'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
		tn.append(_tn)

	run_number = 0
	event_number = 0
	pbar = tqdm.tqdm(range(args.nev))
	for i in pbar:
		if not pythia.next():
			pbar.update(-1)
			continue

		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		jets_hv = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected)))
		if len(jets_hv) < 1:
			pbar.update(-1)
			continue

		event_number = event_number + 1

		jets_h = None
		# do your things with jets here...
		for i, w in enumerate(qweights):
			if w > 0:
				_pthermal = []
				for p in parts_pythia_h_selected:
					parts = thm[i].thermalize(p.perp(), p.eta(), p.phi(), p.m())
					__ = [_pthermal.append(fj.PseudoJet(_p.px(), _p.py(), _p.pz(), _p.e())) for _p in parts] 
				# print(i, 'len', len(_pthermal), 'from', len(parts_pythia_h_selected))
				jets_h = fj.sorted_by_pt(jet_selector(jet_def(_pthermal)))
				for _p in _pthermal:
					tn[i].Fill(run_number, event_number, _p.perp(), _p.eta(), _p.phi(), 111)
			else:
				jets_h = jets_hv

			for j in jets_h:
				hpt[i].Fill(j.perp())
				if j.perp() > 100 and j.perp() < 125:
					for c in j.constituents():
						hz[i].Fill(c.perp() / j.perp())

	rout.cd()
	rout.Write()
	rout.Close()

	pythia.stat()
	pythia.settings.writeFile(args.py_cmnd_out)


if __name__ == '__main__':
	main()
