#!/usr/bin/env python3

import argparse
import os
import pyjetty.alihfjets.hf_data_io as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter
from pyjetty.mputils import JetAnalysisWithRho

import fastjet as fj
import fjcontrib
import fjext

import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)
		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		# self.hinvmass = ROOT.TH1F('hinvmass', 'hinvmass', 400, 1.5, 2.5)
		# self.hinvmass.Sumw2()
		# self.hinvmasspt = ROOT.TH2F('hinvmasspt', 'hinvmasspt', 400, 1.5, 2.5, 50, 2, 12)
		# self.hinvmasspt.Sumw2()
		self.tw = treewriter.RTreeWriter(tree_name='d0', fout=self.fout)
		# jet stuff

		max_eta = 0.9
		self.parts_selector = fj.SelectorPtMin(0.15) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta)

		jet_R0 = 0.4
		self.jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorPtMax(100.0) & fj.SelectorAbsEtaMax(max_eta - 1.05 * jet_R0)
		self.jarho = JetAnalysisWithRho(jet_R=jet_R0, jet_algorithm=fj.antikt_algorithm, particle_eta_max=max_eta)

	def analysis(self, df):
		if len(df) <= 0:
			return

		_parts = self.parts_selector(self.fj_parts)

		self.jarho.analyze_event(_parts)
		_jets = self.jet_selector(self.jarho.jets)

		self.tw.fill_branch('j', _jets)
		self.tw.fill_branch('D0cand', self.fj_Dcands)

		_dRs = [[j.delta_R(dc) for dc in self.fj_Dcands] for j in _jets]
		self.tw.fill_branch('jDdR', _dRs, do_enumerate=True)

		for c in df.columns:
			if 'Particle' in c:
				continue
			self.tw.fill_branch(c, df[c].values)

		self.tw.fill_tree()
				
	def finalize(self):
		self.fout.Write()
		self.fout.Close()
		pinfo(self.fout.GetName(), 'written.')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()

	hfaio = hfdio.HFAnalysisIO()

	hfa = HFAnalysisInvMass(name = args.output)
	hfa.add_selection_range('pt_cand', 2, 1e3)
	hfa.add_selection_range_abs('z_vtx_reco', 10)
	hfa.add_selection_range('pt_prong0', 0.5, 1e3)
	hfa.add_selection_range('pt_prong1', 0.5, 1e3)
	hfa.add_selection_range_abs('eta_cand', 0.8)
	hfa.add_selection_range_abs('nsigTPC_Pi_0', 3)
	hfa.add_selection_range_abs('nsigTPC_Pi_1', 3)
	hfa.add_selection_range_abs('nsigTPC_K_0', 3)
	hfa.add_selection_range_abs('nsigTPC_K_1', 3)
	hfa.add_selection_range_abs('nsigTOF_Pi_0', 3)
	hfa.add_selection_range_abs('nsigTOF_Pi_1', 3)
	hfa.add_selection_range_abs('nsigTOF_K_0', 3)
	hfa.add_selection_range_abs('nsigTOF_K_1', 3)
	# hfa.add_selection_range('d_len', 0.05, 1e3)
	# hfa.add_selection_range('dca', 0.02, 1e3)

	#hfa.add_selection_range('cos_p_xy', 0, 0.99)
	#hfa.add_selection_range('d_len_xy', 0.1, 1e3)
	#hfa.add_selection_range('dca', 0.01, 1e3)

	# hfa.add_selection_range_abs('imp_par_prong0', -0.01)
	# hfa.add_selection_range_abs('imp_par_prong1', -0.01)

	hfaio.add_analysis(hfa)

	# hfaio.load_file("./AnalysisResults.root")
	# hfaio.execute_analyses()
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)

	hfa.finalize()