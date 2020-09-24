#!/usr/bin/env python3

from pyjetty.mputils import pwarning, pinfo, perror
import uproot
import pandas as pd
import numpy as np
import argparse
import os
import tqdm

def pd_tree(path, tname, squery=None):
	try:
		tree = uproot.open(path)[tname]
	except:
		pwarning('error getting', tname, 'from file:', path)
		return None
	if not tree:
		perror('Tree {} not found in file {}'.format(tname, path))
		return None
	df = tree.pandas.df()
	if squery:
		#df.query(squery, inplace=True)
		df = df.query(squery)
		df.reset_index(drop=True)
	return df


def count_D_in_file_merge(fname, _d0cuts_kpi):
	d0_tree_name='PWGHF_TreeCreator/tree_D0'
	event_tree_name='PWGHF_TreeCreator/tree_event_char'
	_ev_cuts = "is_ev_rej == 0 & abs(z_vtx_reco) < 10."
	event_df = pd_tree(path=fname, tname=event_tree_name, squery=_ev_cuts)
	if event_df is None:
		return False
	pinfo('Nev', fname, len(event_df.index))

	d0_df = pd_tree(path=fname, tname=d0_tree_name, squery=_d0cuts_kpi)
	if d0_df is None:
		return False
	pinfo('ND0', fname, len(d0_df.index))
	# pinfo(list(event_df))
	if 'ev_id_ext' in list(event_df):
		d0ev_df = pd.merge(d0_df, event_df, on=['run_number', 'ev_id', 'ev_id_ext'])
	else:
		d0ev_df = pd.merge(d0_df, event_df, on=['run_number', 'ev_id'])
	d0ev_df.query(_ev_cuts, inplace=True)
	d0ev_df_grouped = d0ev_df.groupby(['run_number','ev_id'])
	pinfo('ND0+EvCuts ', fname, len(d0ev_df.index))
	pinfo('GR[ND0+EvCuts] ', fname, len(d0ev_df_grouped))


def count_D_in_file_merge(fname, _d0cuts_kpi):
	d0_tree_name='PWGHF_TreeCreator/tree_D0'
	d0_df = pd_tree(path=fname, tname=d0_tree_name, squery=_d0cuts_kpi)
	if d0_df is None:
		return -1
	d0ev_df_grouped = d0_df.groupby(['run_number','ev_id'])
	# pinfo('ND0 ', fname, len(d0ev_df_grouped))
	return len(d0ev_df_grouped)


def main():
	_d0cuts_base = "(pt_cand > 3.0 & pt_prong0 > 0.6 & pt_prong1 > 0.6 & abs(eta_cand) < 0.8) & "
	_d0cuts_extra = "(dca)<0.03 & abs(cos_t_star)<0.8 & (imp_par_prod)<-0.0001 & (cos_p)>0.9 & "
	_d0cuts_kpi = _d0cuts_base + _d0cuts_extra
	_d0cuts_kpi += "((abs(nsigTPC_Pi_0) < 3. & (abs(nsigTOF_Pi_0) < 3. | nsigTOF_Pi_0 < -900) & abs(nsigTPC_K_1) < 3. & (abs(nsigTOF_K_1) < 3. | nsigTOF_K_1 < -900)) | "
	_d0cuts_kpi += "(abs(nsigTPC_Pi_1) < 3. & (abs(nsigTOF_Pi_1) < 3. | nsigTOF_Pi_1 < -900) & abs(nsigTPC_K_0) < 3. & (abs(nsigTOF_K_0) < 3. | nsigTOF_K_0 < -900)))"

	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='single root file or a file with a list of files to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="prefix output file names", type=str, default='./count_D0.csv')
	args = parser.parse_args()

	fname = args.flist
	if '.root' in fname:
		# count_D_in_file(fname, _d0cuts_kpi)
		nD0s = count_D_in_file_merge(fname, _d0cuts_kpi)
		pinfo(fname, 'N of events with selected D0cand', nD0s)
	else:
		pinfo('reading file list from', fname)
		with open(fname) as f:
			flist = [f.rstrip('\n') for f in f.readlines()]
		pinfo('number of files', len(flist))
		counts = []
		# for ifn, fn in enumerate(flist):
		if args.nfiles > 0:
			flist = flist[:args.nfiles]
		for fn in tqdm.tqdm(flist):
			# pinfo('file', ifn, 'of', len(flist))
			# count_D_in_file(fn, _d0cuts_kpi)
			nD0s = count_D_in_file_merge(fn, _d0cuts_kpi)
			counts.append([fn, nD0s])
		counts_sorted = sorted(counts, key=lambda c: c[1], reverse=True)
		df = pd.DataFrame(counts_sorted, columns=['fname', 'ND0_cand_events'])
		df.to_csv(args.output, index=False)
		pinfo(args.output, 'written.')

if __name__ == '__main__':
	main()
