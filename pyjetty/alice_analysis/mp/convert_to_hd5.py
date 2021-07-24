#!/usr/bin/env python

import cProfile
import uproot
import awkward as ak
import pandas as pd
import argparse

import fastjet as fj
import fjext

import tqdm

class ALICEDataConfig:
	event_tree_name = "PWGHF_TreeCreator/tree_event_char"
	track_tree_name = "PWGHF_TreeCreator/tree_Particle"
	def __init__(self) -> None:
		pass	

# def make_df(df):
# 	d = dict()
# 	# d["particles"] = fjext.vectorize_pt_eta_phi(df['ParticlePt'].values, df['ParticleEta'].values, df['ParticlePhi'].values)
# 	d["parts"] = []
# 	for i in range(len(df['ParticlePt'].values)):
# 		d["parts"].append((df['ParticlePt'].values[i],
# 		             df['ParticleEta'].values[i], 
# 								 df['ParticlePhi'].values[i]))
# 	for dn in df.columns:
# 		if 'Particle' in dn:
# 			continue
# 		d[dn] = df[dn].values[0]
# 	return pd.DataFrame(d)


def make_df_not(df, dfout, args):
	d = dict()
	for dn in df.columns:
		if 'Particle' in dn:
			d[dn] = df[dn].values
		else:
			d[dn] = df[dn].values[0]
	# print(d)
	if dfout is None:
		dfout = pd.DataFrame(d)
		return dfout
	else:
		dfout = dfout.append(pd.DataFrame(d))
	return dfout

def make_df(df, pbar):
	d = dict()
	for dn in df.columns:
		if "run_number" in dn:
			continue
		if "ev_id" in dn:
			continue
		if 'Particle' in dn:
			d[dn] = df[dn].values
			# print(dn, d[dn])
		else:
			d[dn] = df[dn].values[0]
			# print(dn, d[dn])
	dfout = pd.DataFrame(d)
	pbar.update(1)
	return dfout


def process_event(df, args, pbar):
	# fjparts = fjext.vectorize_pt_eta_phi(df['ParticlePt'].values, df['ParticleEta'].values, df['ParticlePhi'].values)
	# jet_def = fj.JetDefinition(fj.antikt_algorithm, 0.4)
	# jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(1.0))
	# particle_selector = fj.SelectorPtMin(0.15) & fj.SelectorAbsRapMax(1.0)
	# fj_particles_selected = particle_selector(fjparts)
	# cs = fj.ClusterSequenceArea(fj_particles_selected, jet_def, jet_area_def)
	# jets = fj.sorted_by_pt(cs.inclusive_jets())
	# if args.debug:
	# 	print('njets: ', len(jets), 'cent: ', df['centrality'].values[0], 'nparts: ', len(fjparts))
	pbar.update(1)

# def analyze_df(df):
# 	for index, row in df:
# 		process_event(DataFrame(row))

def convert(args):
	with uproot.open(args.input)[ALICEDataConfig.event_tree_name] as event_tree:
		event_tree.show()
		print(event_tree.branches)
		#print(event_tree.arrays)
		event_df_orig = event_tree.arrays(library="pd")
		event_df = event_df_orig.query('is_ev_rej == 0')
		event_df.reset_index(drop=True)
		print(event_df)

		with uproot.open(args.input)[ALICEDataConfig.track_tree_name] as track_tree:
			track_tree.show()
			track_df_orig = track_tree.arrays(library="pd")
			# Merge event info into track tree
			track_df = pd.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
			track_df.reset_index(drop=True)
			# (i) Group the track dataframe by event
			#     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
			# # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
			# # df_events = track_df.groupby(['run_number', 'ev_id']).apply(make_df)

			pbar = tqdm.tqdm()
			gby = track_df.groupby(['run_number', 'ev_id'])
			if args.read:
				_df = gby.apply(process_event, args, pbar)
				pbar.close()
			else:
				dfout = None
				# dfout = gby.apply(make_df, dfout, args)
				dfout = gby.apply(make_df, pbar)
				pbar.close()
				dfout.reset_index(drop=True)
				if args.debug:
					print('writing', args.input+'.h5')
				dfout.to_hdf(args.input+'.h5', 'data', mode='a', complevel=9)
				# dfout.to_hdf(args.input+'.h5', 'data', mode='a', complevel=9, format='fixed', data_columns=True)


def read(args):
	df = pd.read_hdf(args.input, 'data')
	pbar = tqdm.tqdm()
	for rn, new_df in df.groupby(["run_number", "ev_id"]):
		process_event(new_df, args, pbar)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', type=str, default='', required=True)
	parser.add_argument('--mc', help="set if this is an MC file",
	                    default=False, action="store_true")
	parser.add_argument('--read', help="read instead of converting", default=False, action="store_true")
	parser.add_argument('--debug', help="flag", default=False, action="store_true")
	parser.add_argument('--profile', help="cProfile", default=False, action="store_true")
	args = parser.parse_args()

	if '.h5' in args.input:
		if args.profile:
			cProfile.run('read(args)')
		else:
			read(args)
	else:
		if args.profile:
			cProfile.run('convert(args)')
		else:
			convert(args)


if __name__ == '__main__':
	main()

