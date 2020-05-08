#!/usr/bin/env python3

import hf_data_io as hfdio
from pyjetty.mputils import perror, pinfo, pwarning
import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)
		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		self.hinvmass = ROOT.TH1F('hinvmass', 'hinvmass', 100, 1, 3)
		self.hinvmass.Sumw2()

	def analysis(self, df):
		for index, row in df.iterrows():
			self.hinvmass.Fill(row['inv_mass'])

	def finalize(self):
		self.fout.Write()
		self.fout.Close()
		pinfo(self.fout.GetName(), 'written.')

if __name__ == '__main__':
	hfaio = hfdio.HFAnalysisIO()

	hfa = HFAnalysisInvMass(name = 'HFAnalysisInvMass')
	hfa.add_selection_range('pt_cand', 2, 1e3)
	hfa.add_selection_range_abs('z_vtx_reco', 10)
	hfa.add_selection_range_abs('eta_cand', 0.8)
	hfa.add_selection_range_abs('nsigTPC_Pi_0', 2)
	hfa.add_selection_range_abs('nsigTPC_Pi_1', 2)
	hfa.add_selection_range_abs('nsigTPC_K_0', 2)
	hfa.add_selection_range_abs('nsigTPC_K_1', 2)
	hfa.add_selection_range_abs('nsigTOF_Pi_0', 2)
	hfa.add_selection_range_abs('nsigTOF_Pi_1', 2)
	hfa.add_selection_range_abs('nsigTOF_K_0', 2)
	hfa.add_selection_range_abs('nsigTOF_K_1', 2)

	hfaio.add_analysis(hfa)

	# hfaio.load_file("./AnalysisResults.root")
	# hfaio.execute_analyses()
	hfaio.execute_analyses_on_file_list('./file_list.txt')

	hfa.finalize()