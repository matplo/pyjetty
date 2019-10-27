#!/usr/bin/env python

from pyjetty.mputils import MPBase, DataIO, JetAnalysis, UniqueString, RTreeWriter
import argparse
import os
import tqdm

import ROOT

import fastjet as fj
import fjcontrib
import fjext

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext


def print_psj(psj, ptcut=0.01):
	if psj.pt() < ptcut:
		return
	area = None
	if psj.has_area():
		area = psj.area()
	print ('pt={} eta={} phi={} usridx={} area={}'.format(psj.pt(), psj.eta(), psj.phi(), psj.user_index(), str(area)))


class GridFastJet(JetAnalysis):
	def __init__(self, **kwargs):
		self.configure_from_args(grid_size=0.1, particle_eta_max=1.)
		self.configure_from_args(**kwargs)
		self._sname = UniqueString.str('h_grid_fast_jet')
		self._nbinsX = int(2 * self.particle_eta_max/self.grid_size)
		self._nbinsY = int(2 * (ROOT.TMath.Pi())/self.grid_size)
		self.grid_eta_phi = ROOT.TH2F(	self._sname, self._sname, 
										self._nbinsX, -self.particle_eta_max, self.particle_eta_max, 
										# int(2 * (ROOT.TMath.Pi())/self.grid_size), -ROOT.TMath.Pi(), ROOT.TMath.Pi())
										self._nbinsY, 0., 2. * ROOT.TMath.Pi())
		super(GridFastJet, self).__init__(**kwargs)

	def parts_from_grid(self):
		for ieta in range(1, self.grid_eta_phi.GetNbinsX()+1):
			for iphi in range(1, self.grid_eta_phi.GetNbinsY()+1):
				if self.grid_eta_phi.GetBinContent(ieta, iphi):
					psj = fj.PseudoJet()
					psj.reset_PtYPhiM(	self.grid_eta_phi.GetBinContent(ieta, iphi), 
										self.grid_eta_phi.GetXaxis().GetBinCenter(ieta),
										self.grid_eta_phi.GetYaxis().GetBinCenter(iphi),
										0.)
					yield psj

	def fill_grid(self, parts):
		self.grid_eta_phi.Reset()
		_tmp = [self.grid_eta_phi.Fill(p.eta(), p.phi(), p.E()) for p in parts]
		self.particles = []
		_tmp = [self.particles.append(p) for p in self.parts_from_grid()]

	def analyze_event(self, parts):
		self.fill_grid(parts)
		JetAnalysis.analyze_event(self, self.particles)


def main(agrs):

	fj.ClusterSequence.print_banner()
	print()

	ja = JetAnalysis(jet_R=args.jet_R, jet_algorithm=fj.antikt_algorithm, particle_eta_max=args.max_eta)
	grfj = GridFastJet(grid_size=0.01, jet_R=args.jet_R, jet_algorithm=fj.antikt_algorithm, particle_eta_max=args.max_eta)
	print(grfj)

	if args.pythia:
		mycfg = []
		pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
		if not pythia:
			print("[e] pythia initialization failed.")
			return

	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()
	tg = RTreeWriter(tree_name='tg', fout=fout)
	ts = RTreeWriter(tree_name='ts', fout=fout)

	for iev in tqdm.tqdm(range(args.nev)):
		if pythia:
			if not pythia.next():
				continue
			parts_pythia = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal, pythiafjext.kCharged])

			ja.analyze_event(parts=parts_pythia)
			if len(ja.jets) < 1:
				continue
			if ja.jets[0].pt() < 10:
				continue
			ts.fill_branch('j', ja.jets)
			ts.fill_branch('p', ja.particles)

			grfj.analyze_event(parts=parts_pythia)
			tg.fill_branch('j', grfj.jets)
			tg.fill_branch('p', grfj.particles)

			ts.fill_tree()
			tg.fill_tree()

	ts.write_and_close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='draw a .draw file', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--pythia', help='run pythia', action='store_true')
	parser.add_argument('--data', help='data events - list of files', type=str, default=None)
	parser.add_argument('--output', default="output.root", type=str)
	parser.add_argument('--jet-R', default=0.4, type=float, help='jet radius')
	parser.add_argument('--max-eta', default=1.0, type=float, help='maximum particle eta')

	args = parser.parse_args()

	main(args)
