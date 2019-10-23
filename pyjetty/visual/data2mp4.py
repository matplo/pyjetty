#!/usr/bin/env python

from __future__ import print_function

import tqdm
import argparse
import os
import numpy as np
import array

import fastjet as fj

from pyjetty.mputils import MPBase
from pyjetty.mputils import DataIO

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT
ROOT.gROOT.SetBatch(True)

class Data2MP4(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt', n_events=500, output='output', run_jetfinder=False)
		super(Data2MP4, self).__init__(**kwargs)
		self.data_io = DataIO(file_list=self.file_list)
		#self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#eta;#varphi', 51, -1, 1, 51, 0., 2.*ROOT.TMath.Pi())
		#self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#eta;#varphi', 51, -1, 1, 51, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		#self.hetaphi_jet = ROOT.TH2F('_hetaphi_jet', '_hetaphi_jet;#eta;#varphi', 51, -1, 1, 51, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		#self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#eta;#varphi', 101, -ROOT.TMath.Pi(), ROOT.TMath.Pi(), 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#eta;#varphi', 51, -1, 1, 51, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		self.hetaphi_jet = ROOT.TH2F('_hetaphi_jet', '_hetaphi_jet;#eta;#varphi', 51, -ROOT.TMath.Pi(), ROOT.TMath.Pi(), 51, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		self.mean_eta = ROOT.TH1F('_hmean_eta', '_hmean_eta;<#eta>', 101, -0.3, 0.3)
		# self.mean_phi = ROOT.TH1F('_hmean_phi', '_hmean_phi;<#varphi>', 101, ROOT.TMath.Pi()-0.3, ROOT.TMath.Pi()+0.3)
		self.mean_phi = ROOT.TH1F('_hmean_phi', '_hmean_phi;<#varphi>', 101, -0.3, +0.3)
		self.mean_e  = ROOT.TH1F('_hmean_e', '_hmean_e;<E (GeV)>', 100, 0, 2)
		self.tc = None
		if self.run_jetfinder:
			fj.ClusterSequence.print_banner()

	def run(self):
		_err_level_tmp = ROOT.gErrorIgnoreLevel
		if self.tc is None:
			self.tc = ROOT.TCanvas('_tc', '_tc', 800, 800)
			self.tc.Divide(2,2)
		ROOT.gErrorIgnoreLevel = ROOT.kWarning
		for iev in tqdm.tqdm(range(self.n_events)):
			parts = self.data_io.load_event(offset=10000)
			if self.run_jetfinder:
				jet_R0 = 0.4
				jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
				jets = fj.sorted_by_pt(jet_def(parts))
				if len(jets) > 0:
					jet = jets[0]
				else:
					jet = fj.PseudoJet()
			else:
				jet = fj.PseudoJet()
			self.hetaphi.Reset()
			self.hetaphi_jet.Reset()
			#_tmp = [self.hetaphi.Fill(p.eta(), jet.delta_phi_to(p), p.e()) for p in parts]
			_tmp = [self.hetaphi.Fill(p.eta(), p.phi()-ROOT.TMath.Pi(), p.e()) for p in parts]
			if self.run_jetfinder:
				# _tmp = [self.hetaphi_jet.Fill(p.eta(), jet.delta_phi_to(p), p.e()) for p in jet.constituents()]
				_tmp = [self.hetaphi_jet.Fill(p.eta(), p.phi()-ROOT.TMath.Pi(), p.e()) for p in jet.constituents()]
			_e = [p.e() for p in parts]
			_phi = [p.phi()-ROOT.TMath.Pi() for p in parts]
			# _phi = [jet.delta_phi_to(p) for p in parts]
			_eta = [p.eta() for p in parts]
			# self.mean_eta.Fill(self.hetaphi.GetMean(1))
			self.mean_eta.Fill(np.mean(_eta))
			# self.mean_phi.Fill(self.hetaphi.GetMean(2))
			self.mean_phi.Fill(np.mean(_phi))
			self.mean_e.Fill(np.mean(_e))
			self.tc.cd()
			# self.hetaphi.Scale(1./self.hetaphi.Integral())
			self.hetaphi.SetMaximum(10.)
			self.hetaphi.SetMinimum(0.)
			self.tc.cd(1)
			self.hetaphi.Draw('colz')
			# self.hetaphi.Draw('lego2')
			self.hetaphi_jet.SetLineColor(ROOT.kRed)
			self.hetaphi_jet.SetFillColor(ROOT.kRed)
			self.hetaphi_jet.SetFillStyle(1001)
			self.hetaphi_jet.Draw('cont3 same')
			self.tc.cd(2)
			self.mean_e.Draw()
			self.tc.cd(3)
			self.mean_phi.Draw()
			self.tc.cd(4)
			self.mean_eta.Draw()
			self.tc.SaveAs('_{}_{}.png'.format(self.output, iev), '.png')
		ROOT.gErrorIgnoreLevel = _err_level_tmp

	def save_mp4(self):
	    # os.system('ffmpeg -i _{}_%01d.png -vcodec mpeg4 -y {}.mp4'.format(self.output, self.output))
	    os.system('ffmpeg -r 4 -f image2 -s 1920x1080 -i _{}_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p {}.mp4'.format(self.output, self.output))

class PythiaJetty(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(pthatmin=100, eta_max=1, jet_pt_min=10, jet_R0=0.4)
		super(PythiaJetty, self).__init__(**kwargs)
		parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=None)
		pyconf.add_standard_pythia_args(parser)
		args = parser.parse_args('')
		args.py_pthatmin = self.pthatmin

		mycfg = []
		self.pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
		if not self.pythia:
			print("[e] pythia initialization failed.")
		self.parts_pythia = None

	def get_event(self):
		if not self.pythia:
			return None
		parts_selector = fj.SelectorAbsEtaMax(self.eta_max)
		jet_selector = fj.SelectorPtMin(self.jet_pt_min) & fj.SelectorAbsEtaMax(self.eta_max - 1.05 * self.jet_R0)
		jet_def = fj.JetDefinition(fj.antikt_algorithm, self.jet_R0)
		while True:
			if not self.pythia.next():
				continue
			self.parts_pythia = pythiafjext.vectorize_select(self.pythia, [pythiafjext.kFinal])
			parts_gen = parts_selector(self.parts_pythia)
			signal_jets = fj.sorted_by_pt(jet_selector(jet_def(self.parts_pythia)))
			if len(signal_jets) < 1:
				continue
			else:
				break
		return self.parts_pythia

class Data2MP4Morph(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt', n_events=500, output='output')
		super(Data2MP4Morph, self).__init__(**kwargs)
		self.data_io = DataIO(file_list=self.file_list)
		# self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#eta;#varphi', 51, -ROOT.TMath.Pi(), ROOT.TMath.Pi(), 51, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		# self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#eta;#varphi', 101, -1, 1, 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi())
		self.hetaphi = ROOT.TH2F('_hetaphi', '_hetaphi;#varphi;#eta', 101, -ROOT.TMath.Pi(), +ROOT.TMath.Pi(), 101, -1, 1)
		self.tc = None
		self.iframe = 0
		self.pythia = PythiaJetty(pthatmin=100, eta_max=1, jet_pt_min=20, jet_R0=0.4)

	def make_a_frame(self, h):
		self.tc.cd()
		for ibx in range(1, h.GetNbinsX() + 1):
			for iby in range(1, h.GetNbinsY() + 1):
				if h.GetBinContent(ibx, iby) < 1e-1:
					h.SetBinContent(ibx, iby, 1e-1)
		h.SetMaximum(100.)
		h.SetMinimum(1e-1)
		# h.Draw('colz')
		h.SetName('{}_f{}'.format(h.GetName().split('_f')[0], self.iframe))
		h.Draw('LEGO2 FB BB')
		
		# h.Draw('surf2 FB BB')
		# h.Draw('CONT3 COLZ')
		# h.Draw('CYL LEGO2 COLZ')
		h.SetLineColor(ROOT.kWhite)
		h.SetLineColorAlpha(ROOT.kWhite, 1.)
		# h.Draw('CYL SURF2')
		# h.Draw('SURF2')
		# h.Draw('PSR LEGO2')
		# self.tc.SetLogz()
		self.tc.SaveAs('_{}_{}.png'.format(self.output, self.iframe), '.png')
		self.iframe = self.iframe + 1

	def morph(self, hfrom, hto, nframes):
		_hmorph = hfrom.Clone("_hmorph")
		_hmorph.SetDirectory(0)
		_hdeltas = hto.Clone('_delta')
		_hdeltas.Add(hfrom, -1.)
		for ifr in tqdm.tqdm(range(1, nframes + 1)):
			fr = 1.0 * ifr / nframes 
			for ibx in range(1, hfrom.GetNbinsX() + 1):
				for iby in range(1, hfrom.GetNbinsY() + 1):
					v = hfrom.GetBinContent(ibx, iby) + fr * _hdeltas.GetBinContent(ibx, iby)
					_hmorph.SetBinContent(ibx, iby, v)
			self.make_a_frame(_hmorph)

	def run(self):
		_err_level_tmp = ROOT.gErrorIgnoreLevel
		if self.tc is None:
			self.tc = ROOT.TCanvas('_tc', '_tc', 800, 800)
		ROOT.gErrorIgnoreLevel = ROOT.kWarning
		_hetaphi_next = self.hetaphi.Clone("_next")
		for iev in tqdm.tqdm(range(self.n_events)):
			parts = self.data_io.load_event(offset=10000)
			_hetaphi_next.Reset()
			#_tmp = [_hetaphi_next.Fill(p.eta(), p.phi()-ROOT.TMath.Pi(), p.e()) for p in parts]
			#_tmp = [_hetaphi_next.Fill(p.phi()-ROOT.TMath.Pi(), p.eta(), p.e()) for p in parts]
			pyev = None
			if self.pythia:
				pyev = self.pythia.get_event()
				if pyev:
					_tmp = [_hetaphi_next.Fill(p.phi()-ROOT.TMath.Pi(), p.eta(), p.e()) for p in pyev]
			if iev > 0:
				self.morph(self.hetaphi, _hetaphi_next, 10)
			self.hetaphi.Reset()
			#_tmp = [self.hetaphi.Fill(p.phi()-ROOT.TMath.Pi(), p.eta(), p.e()) for p in parts]
			if pyev:
				_tmp = [self.hetaphi.Fill(p.phi()-ROOT.TMath.Pi(), p.eta(), p.e()) for p in pyev]
			self.make_a_frame(self.hetaphi)
		ROOT.gErrorIgnoreLevel = _err_level_tmp

	def save_mp4(self):
	    # os.system('ffmpeg -r 10 -i _{}_%d.png -vcodec mpeg4 -y {}.mp4'.format(self.output, self.output))
	    # good page: https://hamelot.io/visualization/using-ffmpeg-to-convert-a-set-of-images-into-a-video/
	    os.system('ffmpeg -r 60 -f image2 -s 1920x1080 -i _{}_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p {}.mp4'.format(self.output, self.output))

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	parser.add_argument('--flist', help='data from a file list', default='', type=str)
	parser.add_argument('--pythia', help='run pythia', default=False, action='store_true')
	parser.add_argument('--output', default="data2mp4", type=str)
	parser.add_argument('--nev', help='number of events', default=500, type=int)
	parser.add_argument('--convert-only', default=False, action='store_true')
	parser.add_argument('--jf', default=False, action='store_true')
	args = parser.parse_args()

	# palette = []
	# paletteSize = 2*2048 #1024
	# for i in range(0, paletteSize):
	# 	palette.append(ROOT.TColor.GetColor(1. - (i * 1.0) / paletteSize, 1. - (i * 1.0) / paletteSize, 1. - (i * 1.0) / paletteSize))
	# palette_i = array.array('i', palette)
	# ROOT.gStyle.SetPalette(paletteSize, palette_i)
	# ROOT.gStyle.SetPalette(paletteSize, palette_i)
	# https://root.cern.ch/doc/master/classTColor.html
	ROOT.gStyle.SetPalette(53)
	# ROOT.gStyle.SetPalette(56)
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)

	if len(args.flist) > 0:
		d2mp4 = Data2MP4(file_list=args.flist, output=args.output, n_events=args.nev, run_jetfinder=args.jf)
	if args.pythia:
		d2mp4 = Data2MP4Morph(file_list=args.flist, output=args.output, n_events=args.nev)
	if d2mp4:
		if args.convert_only is False:
			if not os.path.isfile(args.flist):
				print('[e] input file does not exists')
				return
			d2mp4.run()
		d2mp4.save_mp4()

if __name__ == '__main__':
	main()