#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import sys

import pythia8
import pythiaext
import pythiafjext

from heppy.pythiautils import configuration as pyconf

from pyjetty.mputils import mputils

import ROOT
ROOT.gROOT.SetBatch(1)
ROOT.gSystem.Load('libpyjetty_rutil')


class MyOutputFile(object):
	def __init__(self, foutnamebase, qweight = 0):
		self.foutnamebase = foutnamebase
		self.current_file_number = 0
		self.current_fname = ''
		self.rout = None
		self.hpt = None
		self.hz = None
		self.tn = None
		self.qweight = qweight

	def close_current(self):
		if self.rout:
			self.rout.cd()
			self.rout.Write()
			self.rout.Close()
			self.rout = None
		
	def new_file(self):
		self.close_current()
		if self.rout is None:
			self.current_file_number += 1
			self.current_fname = '{}_w{}_f{}.root'.format(os.path.splitext(self.foutnamebase)[0], self.qweight, self.current_file_number, '.root')
			self.rout = ROOT.TFile(self.current_fname, 'recreate')
			hname = 'hpt_{}'.format(self.qweight)
			htitle = 'hpt w={}'.format(self.qweight)
			self.rout.cd()

			# h = ROOT.TH1F(hname, htitle, 10, mputils.logbins(10, 500, 10))
			self.hpt = ROOT.TH1F(hname, htitle, 10, 0, 250)
			hname = 'hpz_{}'.format(self.qweight)
			htitle = 'hpz w={}'.format(self.qweight)

			self.hz = ROOT.TH1F(hname, htitle, 10, 0, 1)
			hname = 'th_{}'.format(self.qweight)
			htitle = 'thf w={}'.format(self.qweight)

			self.tn = ROOT.TNtuple(	'tree_Particle_gen', 
								'particles from thermalizer {}'.format(self.qweight), 
							   	'run_number:ev_id:ParticlePt:ParticleEta:ParticlePhi:ParticlePID')
		else:
			print('[e] unable to open new file - previous one still non None', file=sys.stderr)
			rout = None
			return False

		if self.rout is None:
			print('[e] unable to open new file - is None', file=sys.stderr)
			return False

		if self.rout.IsOpen() is False:
			print('[e] unable to open new file .IsOpen() is False', file=sys.stderr)
			return False
		return True


class MyOutput(object):
	def __init__(self, foutnamebase, qweights = [0,2,4,6,10]):
		self.qweights = []
		__ = [self.qweights.append(w) for w in qweights]
		_pairs = []
		for w in self.qweights:
			_f = MyOutputFile(foutnamebase, w)
			_pairs.append((w,_f))
		self.files = dict(_pairs)

	def new_files(self):
		rvalue = True
		for w in self.qweights:
			rvalue = rvalue and self.files[w].new_file()
		return rvalue

	def close(self):
		for w in self.qweights:
			self.files[w].close_current()

	def list_files(self):
		_list = []
		for w in self.qweights:
			_list.append(self.files[w].current_fname)
		return _list


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

	qweights = [0, 2, 4, 6, 10] # n-parts quenched

	output = MyOutput('gen_quench_out.root', qweights)
	if output.new_files():
		print('[i] new files with', output.list_files())
	else:
		return

	thm = []
	for i, w in enumerate(qweights):
		_t = ROOT.RUtil.Thermalizer(0.7, w, 1.0, max_eta_hadron)
		thm.append(_t)

	run_number = 0
	event_number = 0
	pbar = tqdm.tqdm(range(args.nev))
	while(pbar.n < args.nev):
		if not pythia.next():
			continue

		parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, False)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

		jets_hv = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected)))
		if len(jets_hv) < 1:
			continue

		if event_number > 0 and event_number % 1000 == 0:
			if output.new_files():
				print('[i] new files:', output.list_files())
			else:
				print('[e] no new file. stop here.')
				break

		event_number = event_number + 1
		pbar.update(1)

		jets_h = None
		# do your things with jets/events here...
		for i, w in enumerate(qweights):
			if w > 0:
				_pthermal = []
				for p in parts_pythia_h_selected:
					parts = thm[i].thermalize(p.perp(), p.eta(), p.phi(), p.m())
					__ = [_pthermal.append(fj.PseudoJet(_p.px(), _p.py(), _p.pz(), _p.e())) for _p in parts] 
				# print(i, 'len', len(_pthermal), 'from', len(parts_pythia_h_selected))
				jets_h = fj.sorted_by_pt(jet_selector(jet_def(_pthermal)))
				for _p in _pthermal:
					output.files[w].tn.Fill(run_number, event_number, _p.perp(), _p.eta(), _p.phi(), 111)
			else:
				jets_h = jets_hv
				for _p in parts_pythia_h_selected:
					output.files[w].tn.Fill(run_number, event_number, _p.perp(), _p.eta(), _p.phi(), 111)

			for j in jets_h:
				output.files[w].hpt.Fill(j.perp())
				if j.perp() > 100 and j.perp() < 125:
					for c in j.constituents():
						output.files[w].hz.Fill(c.perp() / j.perp())

	pbar.close()
	output.close()

	pythia.stat()
	pythia.settings.writeFile(args.py_cmnd_out)


if __name__ == '__main__':
	main()
