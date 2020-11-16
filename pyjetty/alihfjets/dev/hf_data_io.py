from pyjetty.mputils import MPBase, pwarning, pinfo, perror
import random
import uproot
import pandas as pd
import os
import tqdm
# from numba import jit
# import numexpr

class DFSelection(MPBase):
	def __init__(self, **kwargs):
		super(DFSelection, self).__init__(**kwargs)
		self.selection = []
		self.df_selection = None
		self.query_strings = []
		self.query_string = ''

	def add_selection_equal(self, what, val):
		self.selection.append([what, val, None, 0])
		self.query_strings.append('({} == {})'.format(what, val))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_range(self, what, minv, maxv):
		self.selection.append([what, minv, maxv, 1])
		self.query_strings.append('({} > {}) & ({} < {})'.format(what, minv, what, maxv))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_range_abs(self, what, val):
		self.selection.append([what, val, None, 2])
		self.query_strings.append('({} > {}) & ({} < {})'.format(what, -val, what, +val))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_cond(self, what, val1, val2):
		self.selection.append([what, val1, val2, 3])
		self.query_strings.append('((({} > {}) & ({} < {})) | ({} < {}))'.format(what, -val2, what, +val2, what, val1))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_nsig(self, tpcp0, tofp0, tpck1, tofk1, tpcp1, tofp1, tpck0, tofk0, val1, val2):
	#	self.selection.append([tpcp0, tofp0, None, 4])
		self.query_strings.append('((abs({}) < {} & (abs({}) < {} | {} < {}) & abs({}) < {} & (abs({}) < {} | {} < {})) | (abs({}) < {} & (abs({}) < {} | {} < {}) & abs({}) < {} & (abs({}) < {} | {} < {})))'.format(tpcp0, val1, tofp0, val1, tofp0, val2, tpck1, val1, tofk1, val1, tofk1, val2, tpcp1, val1, tofp1, val1, tofp1, val2, tpck0, val1, tofk0, val1, tofk0, val2))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def compile_selection(self, df):
		self.df_selection = True
		for c in self.selection:
			if c[3] == 0:
				self.df_selection = (self.df_selection) & (df[c[0]] == c[1])
			if c[3] == 1:
				self.df_selection = (self.df_selection) & (df[c[0]] > c[1]) & (df[c[0]] < c[2])
			if c[3] == 2:
				self.df_selection = (self.df_selection) & (df[c[0]] > -c[1]) & (df[c[0]] < c[1])


class HFAnalysis(MPBase):
	def __init__(self, **kwargs):
		super(HFAnalysis, self).__init__(**kwargs)
		self.callback = None
		self.event_selection = DFSelection()
		self.event_selection.add_selection_equal('is_ev_rej', 0)
		self.d0_selection = DFSelection()

	def analyze(self, df):
		self.compile_selection(df)
		_df = df[self.df_selection]
		self.analysis(_df)
		if self.callback is not None:
			self.callback(df['ev_id'].values[0])

	def exec_analysis_d0_df(self, df):
		self.pbar.update(1)
		_n_d0s = len(df)
		if _n_d0s < 1:
			return
		if 'ev_id_ext' in list(self.event_df):
			_ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'].values[0], df['ev_id'].values[0], df['ev_id_ext'].values[0])
		else:
			_ev_query = "run_number == {} & ev_id == {}".format(df['run_number'].values[0], df['ev_id'].values[0])
		self.df_tracks = self.track_df.query(_ev_query)
		self.df_tracks.reset_index(drop=True)
		self.analysis(df)
		self.df_tracks = None

	def analyze_slower(self):
		_d0_df = self.d0_df.query(self.d0_selection.query_string)
		pinfo('N d0s ', len(_d0_df.index))
		if 'ev_id_ext' in list(self.event_df):
			d0ev_df = pd.merge(_d0_df, self.event_df, on=['run_number', 'ev_id', 'ev_id_ext'])
		else:
			d0ev_df = pd.merge(_d0_df, self.event_df, on=['run_number', 'ev_id'])
		d0ev_df.query(self.event_selection.query_string, inplace=True)
		d0ev_df_grouped = d0ev_df.groupby(['run_number','ev_id'])
		pinfo('d0s after event cuts ', len(d0ev_df.index))
		pinfo('N d0 groups after event cuts ', len(d0ev_df_grouped))
		with tqdm.tqdm(total=len(d0ev_df_grouped)) as self.pbar:
			_tmp = d0ev_df_grouped.apply(self.exec_analysis_d0_df)
		self.pbar.close()

	# analysis on the single data frame
	# this is something specific to user - overload this one
	def analysis(self, df):
		if len(df) > 0:
			print (df)


class HFAnalysisIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(d0_tree_name='PWGHF_TreeCreator/tree_D0', 
								 track_tree_name='PWGHF_TreeCreator/tree_Particle',
								 event_tree_name='PWGHF_TreeCreator/tree_event_char')
		super(HFAnalysisIO, self).__init__(**kwargs)
		self.analyses = []
		self.d0_df_grouped = None

	def reset_analyses_list(self):
		self.analyses = []

	def add_analysis(self, a):
		self.analyses.append(a)

	def pd_tree(self, path, tname, squery=None):
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

	def load_file(self, path):
		if not os.path.exists(path):
			pwarning('[w] file', path, 'does not exists.')
			return

		self.event_df = self.pd_tree(path, self.event_tree_name)
		if self.event_df is None:
			return False

		self.d0_df = self.pd_tree(path, self.d0_tree_name)
		if self.d0_df is None:
			return False

		self.track_df = self.pd_tree(path, self.track_tree_name)
		if self.track_df is None:
			return False

		return True

	def execute_analyses(self):
		for ana in self.analyses:
			# update_df_references
			ana.event_df 	= self.event_df
			ana.d0_df 		= self.d0_df
			ana.track_df 	= self.track_df
			ana.analyze_slower()
			ana.event_df 	= None
			ana.d0_df 		= None
			ana.track_df 	= None
		# drop the dfs
		self.event_df = None
		self.d0_df = None
		self.track_df = None

	def update_status(self, mark):
		if mark != self.pbar2_mark:
			self.pbar2_mark = mark
			self.pbar2.update(1)

	def execute_analyses_on_file_list(self, file_list, nfiles=0):
		print()
		if os.path.exists(file_list):
			with open(file_list) as f:
				files = f.readlines()
			if int(nfiles) > 0:
				files = files[:nfiles]
			for f in files:
				fn = f.strip('\n')
				pinfo('+file:', fn)
			for f in tqdm.tqdm(files):
				fn = f.strip('\n')
				pinfo('file:', fn)
				if self.load_file(fn):
					self.execute_analyses()
		else:
			perror('file list does not exist', file_list)
		pinfo('done.')

	def __def__(self):
		self.d0_df_grouped = None
