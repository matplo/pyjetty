#!/usr/bin/env python

import argparse
import os
import tqdm
import inspect

from pyjetty.mputils import RTreeReader, find_files, MPBase

import ROOT

class HJetHistogramsBase(MPBase):
	tree_names = ['TT_6_7', 'TT_20_30']
	mean_mV0 = 30.
	def __init__(self, **kwargs):
		self.configure_from_args(itt=-1, pthat=-1, outnamebase='h_jet_histograms')
		super(HJetHistogramsBase, self).__init__(**kwargs)
		self.nevents = 0
		self.init_output()

	def init_output(self):
		pass

	def fill_event(self, tr, itt, pthat):
		pass

	def rescale(self):
		self.fout.cd()
		sigma = self.prof_sigma.GetMean(2)
		print ('[i] {} sigma={} nevents={}'.format(self.fout.GetName(), sigma, self.nevents))

	def write(self):
		self.rescale()
		self.fout.Write()
		print('    written', self.fout.GetName())

class HJetHistogramsEvent(HJetHistogramsBase):
	def __init__(self, **kwargs):
		self.configure_from_args(pthat=-1, outnamebase='h_jet_histograms')
		super(HJetHistogramsEvent, self).__init__(**kwargs)

	def init_output(self):
		self.foutname = self.outnamebase.replace('.root', '') + '_event_pthat_{}.root'.format(self.pthat)
		self.fout = ROOT.TFile(self.foutname, 'recreate')
		# print('[i] opened a new file', self.foutname)
		self.fout.cd()
		self.h_mV0 = ROOT.TH1F('h_mV0', 'h_mV0', 30, 0, 300)
		self.h_mV0_rel = ROOT.TH1F('h_mV0_rel', 'h_mV0_rel', 30, 0, 300 / self.mean_mV0)
		self.prof_sigma = ROOT.TProfile("prof_sigma", "prof_sigma", 2, 0, 1)
		self.prof_sigma_pthat = ROOT.TProfile("prof_sigma_pthat", "prof_sigma_pthat", 1100, 0, 1100)

	def fill_event(self, tr, itt, pthat):
		if self.pthat != pthat:
			return False
		self.nevents = self.nevents + 1
		self.prof_sigma.Fill(0, tr.sigma[0])
		self.prof_sigma_pthat.Fill(self.pthat, tr.sigma[0])
		self.h_mV0.Fill(tr.mV0[0])
		self.h_mV0_rel.Fill(tr.mV0[0]/self.mean_mV0)
		return True

	def rescale(self):
		self.fout.cd()
		sigma = self.prof_sigma.GetMean(2)
		print ('[i] {} sigma={} nevents={}'.format(self.fout.GetName(), sigma, self.nevents))
		hl = [a for a in dir(self) if a.startswith('h_')]
		for hn in hl:
			# print(hn)
			h = getattr(self, hn)
			h.Sumw2()
			hw = h.Clone(h.GetName() + "_weighted")
			hw.Scale(sigma)
			hw.Write()

class HJetHistograms(HJetHistogramsBase):
	def __init__(self, **kwargs):
		self.configure_from_args(itt=-1, pthat=-1, outnamebase='h_jet_histograms')
		super(HJetHistograms, self).__init__(**kwargs)

	def init_output(self):
		self.foutname = self.outnamebase.replace('.root', '') + '{}_pthat_{}.root'.format(self.tree_names[self.itt], self.pthat)
		self.fout = ROOT.TFile(self.foutname, 'recreate')
		# print('[i] opened a new file', self.foutname)
		self.fout.cd()
		self.h_TT = ROOT.TH1F("h_TT", "h_TT;TT p_{T}", 30, 0, 30)
		self.h_TTr = ROOT.TH1F("h_TTr", "h_TTr;TT interval;counts", 2, 0, 2)

		self.h_TT_HM = ROOT.TH1F("h_TT_HM", "h_TT;TT p_{T}", 30, 0, 30)
		self.h_TTr_HM = ROOT.TH1F("h_TTr_HM", "h_TTr;TT interval;counts", 2, 0, 2)

		self.h_TTpThat = ROOT.TH2F("h_TTpThat", "h_TTpThat;TT p_{T};#hat{p_{T}}", 30, 0, 30, 1500, 0, 1500)
		self.h_jetpThat = ROOT.TH2F("h_jetpThat", "h_jetpThat;jet p_{T};#hat{p_{T}}", 30, 0, 60, 1500, 0, 1500)
		self.h_jetpT = ROOT.TH1F("h_jetpT", "h_jetpT;jet p_{T}", 30, 0, 60)
		self.h_jetpT_HM = ROOT.TH1F("h_jetpT_HM", "h_jetpT;jet p_{T}", 30, 0, 60)

		self.h_dphi0 = ROOT.TH1F("h_dphi0", "h_dphi0", 26, 1.5, ROOT.TMath.Pi())
		self.h_dphi1 = ROOT.TH1F("h_dphi1", "h_dphi1", 26, 1.5, ROOT.TMath.Pi())
		self.h_dphi2 = ROOT.TH1F("h_dphi2", "h_dphi2", 26, 1.5, ROOT.TMath.Pi())
		self.h_dphi3 = ROOT.TH1F("h_dphi3", "h_dphi3", 26, 1.5, ROOT.TMath.Pi())

		self.h_dphi0_HM = ROOT.TH1F("h_dphi0_HM", "h_dphi0_HM", 36, 0, ROOT.TMath.Pi())
		self.h_dphi1_HM = ROOT.TH1F("h_dphi1_HM", "h_dphi1_HM", 36, 0, ROOT.TMath.Pi())
		self.h_dphi2_HM = ROOT.TH1F("h_dphi2_HM", "h_dphi2_HM", 36, 0, ROOT.TMath.Pi())
		self.h_dphi3_HM = ROOT.TH1F("h_dphi3_HM", "h_dphi3_HM", 36, 0, ROOT.TMath.Pi())

		self.h_mV0 = ROOT.TH1F('h_mV0', 'h_mV0', 30, 0, 300)
		self.h_mV0_rel = ROOT.TH1F('h_mV0_rel', 'h_mV0_rel', 30, 0, 300 / self.mean_mV0)

		self.prof_sigma = ROOT.TProfile("prof_sigma", "prof_sigma", 2, 0, 1)
		self.prof_sigma_pthat = ROOT.TProfile("prof_sigma_pthat", "prof_sigma_pthat", 1100, 0, 1100)
		# print(self)

	def fill_event(self, tr, itt, pthat):
		if self.itt != itt or self.pthat != pthat:
			return False
		self.nevents = self.nevents + 1
		self.prof_sigma.Fill(0, tr.sigma[0])
		self.prof_sigma_pthat.Fill(self.pthat, tr.sigma[0])
		if max(tr.jet_pt) > self.pthat * 4.:
			return True
		self.h_mV0.Fill(tr.mV0[0])
		self.h_mV0_rel.Fill(tr.mV0[0]/self.mean_mV0)
		self.h_TT.Fill(tr.t_pt[0])
		self.h_TTr.Fill(self.itt)
		self.h_TTpThat.Fill(tr.t_pt[0], tr.pthard[0])

		is_HM = (tr.mV0[0] > 5. * self.mean_mV0 and tr.mV0[0] < 9. * self.mean_mV0)
		if is_HM:
			self.prof_sigma.Fill(1, tr.sigma[0])
			self.h_TT_HM.Fill(tr.t_pt[0])
			self.h_TTr_HM.Fill(self.itt)

		for ij in range(tr.jet_pt.size()):
			if tr.jet_a[ij] < 0.3:
				continue
			self.h_jetpT.Fill(tr.jet_pt[ij])
			self.h_jetpThat.Fill(tr.jet_pt[ij], tr.pthard[0])
			_dphi = abs(tr.t_dphi[ij])
			if tr.jet_pt[ij] > 10 and tr.jet_pt[ij] < 15:
				self.h_dphi0.Fill(_dphi)
			if tr.jet_pt[ij] > 15 and tr.jet_pt[ij] < 20:
				self.h_dphi1.Fill(_dphi)
			if tr.jet_pt[ij] > 20 and tr.jet_pt[ij] < 30:
				self.h_dphi2.Fill(_dphi)
			if tr.jet_pt[ij] > 40 and tr.jet_pt[ij] < 60:
				self.h_dphi3.Fill(_dphi)
			if is_HM:
				self.h_jetpT_HM.Fill(tr.jet_pt[ij])
				if tr.jet_pt[ij] > 10 and tr.jet_pt[ij] < 15:
					self.h_dphi0_HM.Fill(_dphi)
				if tr.jet_pt[ij] > 15 and tr.jet_pt[ij] < 20:
					self.h_dphi1_HM.Fill(_dphi)
				if tr.jet_pt[ij] > 20 and tr.jet_pt[ij] < 30:
					self.h_dphi2_HM.Fill(_dphi)
				if tr.jet_pt[ij] > 40 and tr.jet_pt[ij] < 60:
					self.h_dphi3_HM.Fill(_dphi)
		return True

	def rescale(self):
		self.fout.cd()
		sigma = self.prof_sigma.GetMean(2)
		print ('[i] {} sigma={} nevents={}'.format(self.fout.GetName(), sigma, self.nevents))
		hl = [a for a in dir(self) if a.startswith('h_')]
		for hn in hl:
			# print(hn)
			h = getattr(self, hn)
			h.Sumw2()
			hw = h.Clone(h.GetName() + "_weighted")
			hw.Scale(sigma)
			hw.Write()

class HJetHistogramsContainer(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(outnamebase='h_jet_histograms', event_wise=False)
		super(HJetHistogramsContainer, self).__init__(**kwargs)
		self.hlist = []

	def fill_event(self, tr, itt, pthat):
		filled = False
		for h in self.hlist:
			filled = h.fill_event(tr, itt, pthat)
			if filled:
				break
		if not filled:
			hnew = None
			if self.event_wise:
				hnew = HJetHistogramsEvent(pthat=pthat, itt=itt, outnamebase=self.outnamebase)
			else:
				hnew = HJetHistograms(pthat=pthat, itt=itt, outnamebase=self.outnamebase)
			self.hlist.append(hnew)
			filled = self.hlist[-1].fill_event(tr, itt, pthat)
		if not filled:
			print('[error] combination {} {} not filled!'.format(itt, pthat))
		return filled

	def write(self):
		for h in self.hlist:
			h.write()


def main(args):
	hjcontainer = HJetHistogramsContainer(outnamebase=args.output)
	hjcontainer_event = HJetHistogramsContainer(outnamebase=args.output, event_wise=True)
	file_list = find_files(args.input_dir, 'h_jet_ch_R04_tranges_6-7_20-30_runid_*_pthatmin_*.0_hard.root')
	print("[i] number of files found", len(file_list))
	RTreeReader.print_errors()
	tree_names = ['hjetT_6_7', 'hjetT_20_30', 'evT']
	for ifile in tqdm.tqdm(range(len(file_list))):
		pthat = float(file_list[ifile].split('pthatmin_')[1].split('_hard')[0])
		# print(file_list[ifile], pthat)
		for itt in tqdm.tqdm(range(len(tree_names))):
			tr = None
			container = None
			if tree_names[itt] == 'evT':
				tr = RTreeReader(	tree_name=tree_names[itt], 
									branches = ['pthard', 'sigma', 'mV0'],
									file_name=file_list[ifile])
				container = hjcontainer_event
			else:
				tr = RTreeReader(	tree_name=tree_names[itt], 
									branches = ['t_pt', 'jet_pt', 'jet_a', 'pthard', 't_dphi', 'sigma', 'mV0'],
									file_name=file_list[ifile])
				container = hjcontainer
			nev = tr.tree.GetEntries()
			if args.nev_per_file > 0:
				if args.nev_per_file < nev:
					nev = args.nev_per_file
			for iev in tqdm.tqdm(range(nev)):
				tr.tree.GetEntry(iev)
				container.fill_event(tr, itt, pthat)
			tr.close()
	print()
	print()
	RTreeReader.print_errors()

	hjcontainer.write()
	hjcontainer_event.write()

	# for i in range(tr.tree.GetEntries()):
	# 	tr.tree.GetEntry(i)
	# 	# print (tr.j_pt.size(), tr.ej_pt.size())
	# 	for ie in range(tr.j_pt.size()):
	# 		# print(tr.j_pt[ie])
	# 		pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='analyze hjet events', prog=os.path.basename(__file__))
	parser.add_argument('-d', '--input-dir', help='input directory name - use with default file name', default='.', type=str)
	parser.add_argument('-o', '--output', help='output file name', default='hjet_histograms.root', type=str)
	parser.add_argument('--nev-per-file', help='n events from a file', default=-1, type=int)
	args = parser.parse_args()	
	main(args)