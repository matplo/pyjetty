from pyjetty.mputils import MPBase, pwarning, pinfo, perror
import random
import uproot
import pandas as pd
import fastjet as fj
import fjext
import os
import tqdm


class HFAnalysis(MPBase):
	def __init__(self, **kwargs):
		super(HFAnalysis, self).__init__(**kwargs)
		self.selection = []
		self.df_selection = None

	def add_selection_equal(self, what, val):
		self.selection.append([what, val, None, 0])

	def add_selection_range(self, what, minv, maxv=None):
		self.selection.append([what, minv, maxv, 1])

	def add_selection_range_abs(self, what, val):
		self.selection.append([what, val, None, 2])

	def analyze(self, df):
		self.df_selection = True
		for c in self.selection:
			if c[3] == 0:
				self.df_selection = (self.df_selection) & (df[c[0]] == c[1])
			if c[3] == 1:
				self.df_selection = (self.df_selection) & (df[c[0]] > c[1]) & (df[c[0]] < c[2])
			if c[3] == 2:
				self.df_selection = (self.df_selection) & (df[c[0]] > -c[1]) & (df[c[0]] < c[1])
		_df = df.loc[self.df_selection]
		self.analysis(_df)

	# analysis on the single data frame
	# this is something specific to user - overload this one
	def analysis(self, df):
		if len(df) > 0:
			print (df)


class HFAnalysisIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(tree_name='PWGHF_TreeCreator/tree_D0')
		self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
		super(HFAnalysisIO, self).__init__(**kwargs)
		self.analyses = []
		self.track_df_grouped = None

	def reset_analyses_list(self):
		self.analyses = []

	def add_analysis(self, a):
		self.analyses.append(a)

	def load_file(self, path):
		if not os.path.exists(path):
			pwarning('[w] file', path, 'does not exists.')
			return
		event_tree = uproot.open(path)[self.event_tree_name]
		if not event_tree:
			print('[e] Tree {} not found in file {}'.format(self.event_tree_name, path))
			return False
		event_df_orig = event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		event_df_orig.reset_index(drop=True)
		event_df = event_df_orig.query('is_ev_rej == 0')
		event_df.reset_index(drop=True)
		# Load track tree into dataframe
		track_tree = uproot.open(path)[self.tree_name]
		if not track_tree:
			print('[e] Tree {} not found in file {}'.format(tree_name, path))
			return False
		track_df_orig = track_tree.pandas.df()
		# Merge event info into track tree
		track_df = pd.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
		self.track_df_grouped = track_df.groupby(['run_number','ev_id'])
		# (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles

	def execute_analyses(self):
		[self.track_df_grouped.apply(a.analyze) for a in self.analyses]

	def execute_analyses_on_file_list(self, file_list):
		if os.path.exists(file_list):
			with open(file_list) as f:
				files = f.readlines()
			for f in tqdm.tqdm(files):
				fn = f.strip('\n')
				pinfo('processing file:', fn)
				self.load_file(fn)
				self.execute_analyses()
		else:
			perror('file list does not exist', file_list)
		pinfo('done.')

	def __def__(self):
		self.track_df_grouped = None
