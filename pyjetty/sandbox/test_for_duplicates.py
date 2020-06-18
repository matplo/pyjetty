#!/usr/bin/env python

from __future__ import print_function

import tqdm
import argparse
import os
import numpy as np
import pandas as pd
import uproot

from pyjetty.mputils import pinfo, pwarning, perror, pindent

def main():
	parser = argparse.ArgumentParser(description='test duplicate entries', prog=os.path.basename(__file__))
	parser.add_argument('fname', help='input file', default='', type=str)
	args = parser.parse_args()

	event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
	event_tree = uproot.open(args.fname)[event_tree_name]
	if not event_tree:
		perror('Tree {} not found in file {}'.format(event_tree_name, args.fname))
		return False

	pinfo(args.fname)
	event_df_orig = event_tree.pandas.df()
	len_event_df_orig = len(event_df_orig)

	df_event_accepted = event_df_orig.query('is_ev_rej == 0')
	df_event_accepted.reset_index(drop=True)
	len_event_df_accepted = len(df_event_accepted)

	event_df_nodup = df_event_accepted.drop_duplicates()
	len_event_df_nodup = len(event_df_nodup)

	if len_event_df_accepted != len_event_df_nodup:
		perror('original event length:', len_event_df_orig, 'accepted:', len_event_df_accepted, 'nodup:', len_event_df_nodup)
	else:
		pindent('original event length:', len_event_df_orig, 'accepted:', len_event_df_accepted, 'nodup:', len_event_df_nodup)

	track_tree_name = 'PWGHF_TreeCreator/tree_Particle'
	track_tree = uproot.open(args.fname)[track_tree_name]
	if not track_tree:
		perror('Tree {} not found in file {}'.format(tree_name, args.fname))
		return False
	track_df_orig = track_tree.pandas.df()
	track_df = pd.merge(track_df_orig, event_df_nodup, on=['run_number', 'ev_id'])
	len_track_df = len(track_df)
	track_df_nodup = track_df.drop_duplicates()
	len_track_df_nodup = len(track_df_nodup)
	if len_track_df_nodup < len_track_df:
		perror('track+event rows:', len_track_df, 'nodup:', len_track_df_nodup)
	else:
		pindent('track+event rows:', len_track_df, 'nodup:', len_track_df_nodup)		
	track_df_grouped = track_df.groupby(['run_number','ev_id'])
	len_track_df_grouped = len(track_df_grouped)
	if len_track_df_grouped <= len_event_df_nodup:
		pindent ('track+event length grouped:', len_track_df_grouped)
	else:
		perror ('track+event length grouped:', len_track_df_grouped)
	# track_df_nodup = track_df_grouped.drop_duplicates()
	# print ('track+event length no dup:', len(track_df_nodup))

	# from James
	# Check if there are duplicated tracks in an event.
	duplicate_selection = ['run_number', 'ev_id', 'ParticlePt', 'ParticleEta', 'ParticlePhi']
	#if use_ev_id_ext:
	#  duplicate_selection.append('ev_id_ext')
	duplicate_rows_df = track_df.duplicated(duplicate_selection)
	n_duplicates = sum(duplicate_rows_df)
	pindent('2nd pass: using duplicate selection ', duplicate_selection)
	if n_duplicates > 0:
		perror('2nd pass: there appear to be {} duplicate particles in the dataframe'.format(n_duplicates))
		perror('this is: {:.2} of all tracks'.format(n_duplicates/len_track_df))
		track_df_nodup = track_df.drop_duplicates(duplicate_selection, inplace=False)
		pwarning('new count rows for particles:', len(track_df_nodup), 'old count:', len_track_df)
	else:
		pindent('no duplicate particles found')

if __name__ == '__main__':
	main()
