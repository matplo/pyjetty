from pyjetty.mputils import MPBase
import random
import uproot
import pandas as pd
import fjext


class DataEvent(object):
	def __init__(self, particles, run_number, ev_id):
		self.particles = particles
		self.run_number = run_number
		self.ev_id = ev_id


class DataFileIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_input = None, tree_name='tree_Particle')
		self.event_number = 0
		self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
		super(DataFileIO, self).__init__(**kwargs)
		self.reset_dfs()
		if self.file_input:
			self.load_file(self.file_input, self.tree_name)

	def reset_dfs(self):
		self.track_df_orig = None
		self.track_df = None
		self.track_df_grouped = None
		self.df_events = None
		self.event_tree = None
		self.event_df_orig = None
		self.event_df = None
		self.track_tree = None
		self.track_tree_name = None

	def get_event(self, df_tracks):
		# Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
		event = DataEvent([], -1, -1)
		event.particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)  
		if len(event.particles) > 0:
			event.run_number = float(df_tracks['run_number'].values[0])
			event.ev_id = float(df_tracks['ev_id'].values[0])
		else:
			event.run_number = -1
			event.ev_id = -1
		return event

	def load_file(self, file_input, tree_name='tree_Particle'):
		self.file_input = file_input
		self.tree_name = tree_name
		self.reset_dfs()
		# Load event tree into dataframe, and apply event selection
		self.event_tree = uproot.open(file_input)[self.event_tree_name]
		if not self.event_tree:
			print('[e] Tree {} not found in file {}'.format(self.event_tree_name, file_input))
			return False
		self.event_df_orig = self.event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		self.event_df_orig.reset_index(drop=True)
		self.event_df = self.event_df_orig.query('is_ev_rej == 0')
		self.event_df.reset_index(drop=True)
		# Load track tree into dataframe
		self.track_tree_name = 'PWGHF_TreeCreator/{}'.format(tree_name)
		self.track_tree = uproot.open(file_input)[self.track_tree_name]
		if not self.track_tree:
			print('[e] Tree {} not found in file {}'.format(tree_name, file_input))
			return False
		self.track_df_orig = self.track_tree.pandas.df()
		# Merge event info into track tree
		self.track_df = pd.merge(self.track_df_orig, self.event_df, on=['run_number', 'ev_id'])
		# (i) Group the track dataframe by event
		#     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
		self.track_df_grouped = self.track_df.groupby(['run_number','ev_id'])
		# (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
		self.df_events = self.track_df_grouped.apply(self.get_event)
		self.event_number = self.event_number + len(self.df_events)
		return True


#random order of files; files do not repeat
class DataIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt')	
		super(DataIO, self).__init__(**kwargs)
		self.current_event_in_file = 0
		self.file_io = None
		self.list_of_files = []
		self.read_file_list()

	def set_file_list(self, sfile):
		self.file_list = sfile
		self.read_file_list()

	def read_file_list(self):
		self.list_of_files = []
		with open(self.file_list) as f:
			self.list_of_files = [fn.strip() for fn in f.readlines()]

	def open_file(self):
		self.file_io = None
		self.current_event_in_file = 0
		if len(self.list_of_files) > 0:
			afile = random.choice(self.list_of_files)
			self.list_of_files.remove(afile)
		else:
			print('[w] no more files to open.')
			return
		print('[i] opening data file', afile)
		self.file_io = DataFileIO(file_input=afile)
		print('    number of events', self.current_file_number_of_events())

	def current_file_number_of_events(self):
		nev = 0
		if self.file_io:
			if self.file_io.df_events is not None:
				nev = len(self.file_io.df_events)
		return nev

	def load_event(self):
		self.particles = None
		if self.file_io is None:
			self.open_file()
		if self.file_io is None:
			print('[e] unable to load the data file')
			return None
		if self.current_event_in_file >= self.current_file_number_of_events():
			self.current_event_in_file = 0
			self.file_io = None
			return self.load_event()
		event = self.file_io.df_events[self.current_event_in_file]
		_tmp = [p.set_user_index(10000+ip) for ip,p in enumerate(event.particles)]
		# print('loaded event:', self.current_event_in_file)
		self.current_event_in_file = self.current_event_in_file + 1
		self.particles = event.particles
		return self.particles


#random order of files; files can repeat
class DataBackgroundIO(DataIO):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt')	
		super(DataBackgroundIO, self).__init__(**kwargs)

	def open_file(self):
		self.file_io = None
		self.current_event_in_file = 0
		afile = random.choice(self.list_of_files)
		print('[i] opening data file', afile)
		self.file_io = DataFileIO(file_input=afile)
		print('    number of events', self.current_file_number_of_events())
