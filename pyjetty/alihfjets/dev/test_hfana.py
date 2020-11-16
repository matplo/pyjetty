#!/usr/bin/env python3

import argparse
import os
import pyjetty.alihfjets.dev.hf_data_io as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter, jet_analysis

import fastjet as fj
import fjext
import fjtools

import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)
		self.tw = treewriter.RTreeWriter(tree_name='d0', file_name = self.name+'.root')
		self.twjc = treewriter.RTreeWriter(name = 'd0jc', file_name = self.name+'_c.root')
		
	def analysis(self, df):

		djmm = fjtools.DJetMatchMaker()
		djmm.set_ch_pt_eta_phi(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values)
		djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
		djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
		djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)

		self.tw.fill_branches(dpsj = djmm.Ds)
		self.tw.fill_tree()

		for id0, d0 in enumerate(djmm.Ds):
			_parts_and_ds = djmm.match(0.005, id0)
			_parts_and_ds.push_back(d0)
			ja = jet_analysis.JetAnalysis(jet_R = 0.4, particle_eta_max=0.9, jet_pt_min=2.0)
			ja.analyze_event(_parts_and_ds)
			if len(ja.jets) < 1:
				continue
			jets = ja.jets_as_psj_vector()
			djets = djmm.filter_D0_jets(jets)
			if len(djets) > 0:
				j = djets[0]
				dcand = djmm.get_Dcand_in_jet(j)
				sja = jet_analysis.JetAnalysis(jet_R = 0.1, particle_eta_max=0.9, jet_pt_min=2.0)
				sja.analyze_event(j.constituents())
				lsj = fj.sorted_by_pt(sja.jets_as_psj_vector())
				if len(lsj) < 1:
					continue
				sj_dcand = djmm.get_Dcand_in_jet(lsj[0])
				is_Dsj = 0
				if len(sj_dcand) > 0:
					# if sj_dcand[0].m() == dcand[0].m() and sj_dcand[0].perp() == dcand[0].perp():
					if sj_dcand[0].delta_R(dcand[0]) == 0.0:
						is_Dsj = 1
				self.twjc.fill_branches(jet = j, dR = j.delta_R(dcand[0]), D = dcand[0], lsj = lsj[0], Dsj = is_Dsj
										, a10 = fjext.angularity(j,  1.0, 0.4)
										, a15 = fjext.angularity(j,  0.5, 0.4)
										, a20 = fjext.angularity(j,  0.0, 0.4)
										, a30 = fjext.angularity(j, -1.0, 0.4))
				self.twjc.fill_tree()

			if len(djets) > 1:
				perror("more than one jet per D candidate?")

		return True

		## for index, row in df.iterrows():
		## 	# self.hinvmass.Fill(row['inv_mass'])
		## 	# self.hinvmasspt.Fill(row['inv_mass'], row['pt_cand'])
		## 	for c in df.columns:
		## 		self.tw.fill_branch(c, row[c])
		## 	self.tw.fill_tree()
				
	def finalize(self):
		self.tw.write_and_close()
		self.twjc.write_and_close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()

	hfaio = hfdio.HFAnalysisIO()

	hfa = HFAnalysisInvMass(name = args.output)
	#hfa.add_selection_range('pt_cand', 2, 1e3)
	hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
	#hfa.add_selection_range('pt_prong0', 0.5, 1e3)
	#hfa.add_selection_range('pt_prong1', 0.5, 1e3)
	hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
	hfa.d0_selection.add_selection_range('dca', -1, 0.03)
	hfa.d0_selection.add_selection_range_abs('cos_t_star', 0.8)
	hfa.d0_selection.add_selection_range('imp_par_prod', -1, -0.0001)
	hfa.d0_selection.add_selection_range('cos_p', 0.9, 3)
	hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)

	hfaio.add_analysis(hfa)

	# hfaio.load_file("./AnalysisResults.root")
	# hfaio.execute_analyses()
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)

	hfa.finalize()
