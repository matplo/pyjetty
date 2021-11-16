#!/usr/bin/env python3

import numpy as np
import argparse
import os
import pyjetty.alihfjets.dev.hfjet.process_io_mc_hf as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter, jet_analysis


import fastjet as fj
import fjext
import fjtools
import fjcontrib

import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)

		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		self.twjc = treewriter.RTreeWriter(tree_name='d0jc', fout=self.fout)
		self.twjc_gen = treewriter.RTreeWriter(tree_name='d0jc_gen', fout=self.fout)
	
	
	def analysis(self, df):
		#print(df)

		m_array = np.full((self.df_tracks['ParticlePt'].values.size), 0.1396)	

		djmm = fjtools.DJetMatchMaker()
		djmm.set_ch_pt_eta_phi_m(self.df_tracks['ParticlePt'].values, self.df_tracks['ParticleEta'].values, self.df_tracks['ParticlePhi'].values, m_array)
		djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
		djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
		djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)


		#run for each D candidate to build jet
		for id0, d0 in enumerate(djmm.Ds):

			#daughter tracks matching
			_parts_and_ds = djmm.match(0.005, id0)
			#replacing daughter tracks with matched D0 candidate
			#including D0 	
			_parts_and_ds.push_back(d0)
	
			#jet reconstruction with D0 and charged particle
			jetR=0.4
			ja = jet_analysis.JetAnalysis(jet_R = jetR,jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9, jet_pt_min=5.0)
			ja.analyze_event(_parts_and_ds)
			if len(ja.jets) < 1:
				continue
			jets = ja.jets_as_psj_vector()

			#filtering D0 jets
			djets = djmm.filter_D0_jets(jets)


			if len(djets) > 0:
				j = djets[0]
				dcand = djmm.get_Dcand_in_jet(j)
				
				#number of constitutents > 1
				#if len(j.constituents())<=1:
				#	continue

				#jets with the winner take all axis################################		
				jet_def_wta = fj.JetDefinition(fj.cambridge_algorithm, 2*jetR)
				jet_def_wta.set_recombination_scheme(fj.WTA_pt_scheme)
				#print('WTA jet definition is:', jet_def_wta)
				reclusterer_wta =  fjcontrib.Recluster(jet_def_wta)
				jet_wta = reclusterer_wta.result(j)
				################################
					
				D_cosp=df['cos_p'].values
				D_cosTStar=df['cos_t_star'].values
				D_NormalisedDecayLength=df['norm_dl_xy'].values	
				D_ImpactParameterProduct=df['imp_par_prod'].values
				D_cand_type = df['cand_type'].values
				
				self.twjc.fill_branches(jet = j,jetWta =jet_wta, dR = j.delta_R(dcand[0]), dRWTA = jet_wta.delta_R(dcand[0]),  D = dcand[0], cos_p = float(D_cosp[id0]), D_cos_t_star = float(D_cosTStar[id0]), D_norm_dlxy=float(D_NormalisedDecayLength[id0]),D_imp_par_prod=float(D_ImpactParameterProduct[id0]),Dmeson_cand_type=float(D_cand_type[id0]))
				self.twjc.fill_tree()
		
			if len(djets) > 1:
				perror("more than one jet per D candidate?")

		return True
	
	def analysis_gen(self, df):

		m_gen_array = np.full((self.df_gen_tracks['ParticlePt'].values.size), 0.1396)

		djmm_gen = fjtools.DJetMatchMaker()
		djmm_gen.set_ch_pt_eta_phi_m(self.df_gen_tracks['ParticlePt'].values, self.df_gen_tracks['ParticleEta'].values, self.df_gen_tracks['ParticlePhi'].values, m_gen_array)
		djmm_gen.set_Ds_pt_eta_phi(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values)

		for id0, d0 in enumerate(djmm_gen.Ds):
			_parts_and_ds = djmm_gen.ch
			_parts_and_ds.push_back(d0)
			ja_gen = jet_analysis.JetAnalysis(jet_R = 0.4,jet_RecombinationScheme=fj.E_scheme, particle_eta_max=0.9, jet_pt_min=5.0)
			ja_gen.analyze_event(_parts_and_ds)
			if len(ja_gen.jets) < 1:
				continue

			jets = ja_gen.jets_as_psj_vector()
			djets = djmm_gen.filter_D0_jets(jets)

			if len(djets) > 0:
				j = djets[0]
				dcand = djmm_gen.get_Dcand_in_jet(j)
				D_cand_type = df['cand_type'].values
				self.twjc_gen.fill_branches(jet = j,dR = j.delta_R(dcand[0]),D = dcand[0],Dmeson_cand_type=float(D_cand_type[id0]))
				self.twjc_gen.fill_tree()

			if len(djets) > 1:
                       		perror("more than one jet per D candidate?")

		return True

	def finalize(self):
	
		self.fout.Write()
		self.fout.Close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()
	
	hfaio = hfdio.HFAnalysisIO()
	
	hfa = HFAnalysisInvMass(name = args.output)

	hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
	hfa.event_selection.add_selection_equal('is_ev_rej', 0)	
	
	hfa.d0_selection.add_selection_range('inv_mass', 1.46, 2.26)
	hfa.d0_selection.add_selection_range('dca', -1, 0.03)
	hfa.d0_selection.add_selection_range_abs('cos_t_star', 1)
	hfa.d0_selection.add_selection_range('pt_prong0', 0.7, 1e3)
	hfa.d0_selection.add_selection_range('pt_prong1', 0.7, 1e3)
	hfa.d0_selection.add_selection_range_abs('imp_par_prong0', 0.1)
	hfa.d0_selection.add_selection_range_abs('imp_par_prong1', 0.1)
	hfa.d0_selection.add_selection_range('imp_par_prod', -1,0.001)
	hfa.d0_selection.add_selection_range('cos_p', 0.8, 1)
	hfa.d0_selection.add_selection_range('cos_p_xy', 0, 1)
	hfa.d0_selection.add_selection_range('norm_dl_xy', 0, 1e3)
	#topomatic cut suggested by D2H.
	#hfa.d0_selection.add_selection_range_abs('max_norm_d0d0exp',2)

	hfa.d0_selection.add_selection_range('pt_cand', 2, 1e3)
	hfa.d0_selection.add_selection_range_abs('eta_cand', 0.8)
	hfa.d0_selection.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)
	

	

	hfaio.add_analysis(hfa)
	
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)
	
	hfa.finalize()
