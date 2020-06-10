from pyjetty.mputils import MPBase, pwarning, pinfo, perror
import random
import uproot
import pandas as pd
import fastjet as fj
import fjext
import os
import tqdm
# from numba import jit
# import numexpr

class HFAnalysis(MPBase):
	def __init__(self, **kwargs):
		super(HFAnalysis, self).__init__(**kwargs)
		self.selection = []
		self.df_selection = None
		self.query_strings = []
		self.query_string = ''
		self.callback = None
		self.fj_parts = None
		self.fj_Dcands = None

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

	def compile_selection(self, df):
		self.df_selection = True
		for c in self.selection:
			if c[3] == 0:
				self.df_selection = (self.df_selection) & (df[c[0]] == c[1])
			if c[3] == 1:
				self.df_selection = (self.df_selection) & (df[c[0]] > c[1]) & (df[c[0]] < c[2])
			if c[3] == 2:
				self.df_selection = (self.df_selection) & (df[c[0]] > -c[1]) & (df[c[0]] < c[1])

	def analyze(self, df):
		self.compile_selection(df)
		_df = df[self.df_selection]
		self.fj_parts = fjext.vectorize_pt_eta_phi(_df['ParticlePt'].values, _df['ParticleEta'].values, _df['ParticlePhi'].values)
		self.fj_Dcands = fjext.vectorize_pt_eta_phi(_df['pt_cand'].values, _df['eta_cand'].values, _df['phi_cand'].values)
		_df = _df.drop(columns=['ParticlePt', 'ParticleEta', 'ParticlePhi'])
		_df.reset_index(drop=True)
		_df.drop_duplicates(inplace=True)
		self.analysis(_df)
		if self.callback is not None:
			self.callback(df['ev_id'].values[0])

	def analyze_slower(self, df):
		_df = df.query(self.query_string)
		self.analysis(_df)

	# analysis on the single data frame
	# this is something specific to user - overload this one
	def analysis(self, df):
		if len(df) > 0:
			print (df)


class DataEvent(object):
	def __init__(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)
	def configure_from_args(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def __str__(self):
		s = []
		s.append('[i] DataEvent')
		for a in self.__dict__:
			sval = str(getattr(self, a))
			if len(sval) > 200:
				sval = sval[:200]
			s.append('   {} = {}'.format(str(a), sval))
		return '\n'.join(s)


class HFAnalysisIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(d0_tree_name='PWGHF_TreeCreator/tree_D0', 
								 track_tree_name='PWGHF_TreeCreator/tree_Particle',
								 event_tree_name='PWGHF_TreeCreator/tree_event_char',
								 enable_jet=True,
								 enable_d0=True,
								 offset_parts=0)
		super(HFAnalysisIO, self).__init__(**kwargs)
		self.analyses = []
		self.df_grouped = None
		self.df_events = None

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
			df = df.query(squery)
			df.reset_index(drop=True)
		return df

	def load_file(self, path):
		self.df_grouped = None
		self.df_events = None
		event_df = self.pd_tree(path=path, tname=self.event_tree_name, squery='is_ev_rej == 0')
		if event_df is None:
			return False
		d0_df = None
		if self.enable_d0:
			d0_df = self.pd_tree(path=path, tname=self.d0_tree_name)
		track_df = None
		track_df_fj = None
		if self.enable_jet:
			track_df = self.pd_tree(path=path, tname=self.track_tree_name)
		df_merged = None
		if d0_df is None:
			df_merged = event_df
		else:
			df_merged = pd.merge(event_df, d0_df, on=['run_number', 'ev_id'])
		if track_df is None:
			pass
		else:
			df_merged = pd.merge(df_merged, track_df, on=['run_number', 'ev_id'])
		self.df_grouped = df_merged.groupby(['run_number','ev_id'])
		# self.df_grouped.describe()
		# self.df_events = self.df_grouped.apply(self.get_event)
		return True

	def get_event(self, df):
		# Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
		_parts = fjext.vectorize_pt_eta_phi(df['ParticlePt'].values, df['ParticleEta'].values, df['ParticlePhi'].values)  
		#event = DataEvent(particles=_parts)
		#for c in df_tracks.columns:
		#	if 'Particle' not in c:
		#		event.__setattr__(c, df_tracks[c].values[0])
		#if len(event.particles) <= 0:
		#	# backward compat
		#	event.run_number = -1
		#	event.ev_id = -1
		_event = pd.DataFrame({"run_number": df["run_number"], "ev_id": df["ev_id"], "ParticleFJ": _parts})
		event = pd.merge(df, _event, on=['run_number', 'ev_id'])
		return event

	def execute_analyses(self):
		[self.df_grouped.apply(a.analyze) for a in self.analyses]
		# [self.df_events.apply(a.analyze) for a in self.analyses]

	def update_status(self, mark):
		if mark != self.pbar2_mark:
			self.pbar2_mark = mark
			self.pbar2.update(1)

	def execute_analyses_on_file_list(self, file_list, nfiles=0):
		self.pbar2 = tqdm.tqdm(mininterval=5, maxinterval=20)
		self.pbar2_mark = None
		for a in self.analyses:
			a.callback = self.update_status
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
				if self.load_file(fn):
					self.execute_analyses()
			self.pbar2.close()
		else:
			perror('file list does not exist', file_list)
		pinfo('done.')

	def __def__(self):
		self.df_grouped = None

# ------------------------------

class HFAnalysisIOnoJet(MPBase):
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
		try:
			event_tree = uproot.open(path)[self.event_tree_name]
		except:
			pwarning('error getting', self.event_tree_name, 'from file:', path)
			return False
		if not event_tree:
			perror('Tree {} not found in file {}'.format(self.event_tree_name, path))
			return False
		event_df_orig = event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		event_df_orig.reset_index(drop=True)
		event_df = event_df_orig.query('is_ev_rej == 0')
		event_df.reset_index(drop=True)
		# Load track tree into dataframe
		try:
			track_tree = uproot.open(path)[self.tree_name]
		except:
			pwarning('error getting', self.tree_name, 'from file:', path)
			return False
		if not track_tree:
			perror('Tree {} not found in file {}'.format(tree_name, path))
			return False
		track_df_orig = track_tree.pandas.df()
		# Merge event info into track tree
		track_df = pd.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
		self.track_df_grouped = track_df.groupby(['run_number','ev_id'])
		return True

	def execute_analyses(self):
		[self.track_df_grouped.apply(a.analyze) for a in self.analyses]

	def update_status(self, mark):
		if mark != self.pbar2_mark:
			self.pbar2_mark = mark
			self.pbar2.update(1)

	def execute_analyses_on_file_list(self, file_list, nfiles=0):
		self.pbar2 = tqdm.tqdm(mininterval=20, maxinterval=60)
		self.pbar2_mark = None
		for a in self.analyses:
			a.callback = self.update_status
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
				if self.load_file(fn):
					self.execute_analyses()
			self.pbar2.close()
		else:
			perror('file list does not exist', file_list)
		pinfo('done.')

	def __def__(self):
		self.track_df_grouped = None
