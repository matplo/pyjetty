from pyjetty.mputils import MPBase
import random
import uproot
import pandas as pd
import fastjet as fj
import fjext


#class DataEvent(MPBase):
#	def __init__(self, **kwargs):
#		kwargs['name'] = 'noUniqueName'
#		super(DataEvent, self).__init__(**kwargs)
#		# self.particles = particles
#		# self.run_number = run_number
#		# self.ev_id = ev_id

# with less overhead
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
		# event_df = self.event_df.loc[(self.event_df['run_number'] == self.run_number) & (self.event_df['ev_id'] == self.ev_id)]
		# event = DataEvent([], -1, -1, self.event_df)
		# event = DataEvent(particles=[], run_number=-1, ev_id=-1)
		# Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
		_parts = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)  
		event = DataEvent(particles=_parts)
		for c in df_tracks.columns:
			if 'Particle' not in c:
				event.__setattr__(c, df_tracks[c].values[0])
		if len(event.particles) <= 0:
			# backward compat
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
		# self.event_df_orig = self.event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		self.event_df_orig = self.event_tree.pandas.df()
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


#random order of files; files do not repeat; load_event - single event return
class DataIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt', tree_name='tree_Particle', random_file_order=True)
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
			if self.random_file_order:
				afile = random.choice(self.list_of_files)
			else:
				afile = self.list_of_files[0]
			self.list_of_files.remove(afile)
		else:
			print('[w] no more files to open.')
			return
		print('[i] opening data file', afile)
		self.file_io = DataFileIO(file_input=afile, tree_name=self.tree_name)
		print('    number of events', self.current_file_number_of_events(), 'in tree', self.tree_name)
		print('    files to go', len(self.list_of_files))

	def current_file_number_of_events(self):
		nev = 0
		if self.file_io:
			if self.file_io.df_events is not None:
				nev = len(self.file_io.df_events)
		return nev

	def load_event(self, offset = 0):
		self.particles = None
		if self.file_io is None:
			self.open_file()
		if self.file_io is None:
			print('[e] unable to load the data file')
			return None
		if self.current_event_in_file >= self.current_file_number_of_events():
			self.current_event_in_file = 0
			self.file_io = None
			return self.load_event(offset=offset)
		self.event = self.file_io.df_events[self.current_event_in_file]
		# print ('reset indexes')
		# _tmp = [p.set_user_index(0) for ip,p in enumerate(event.particles)]
		# print('loaded event:', self.current_event_in_file)
		self.current_event_in_file = self.current_event_in_file + 1
		self.particles = fj.vectorPJ()
		for ip,p in enumerate(self.event.particles):
			p.set_user_index(offset + ip)
			self.particles.push_back(p)
		return self.particles

	def open_afile(self, afile):
		if self.file_io:
			if self.file_io.file_input == afile:
				return
		self.file_io = None
		self.current_event_in_file = 0
		print('[i] opening (a) data file', afile)
		self.file_io = DataFileIO(file_input=afile, tree_name=self.tree_name)
		print('    number of events', self.current_file_number_of_events(), 'in tree', self.tree_name)

	def load_event_with_loc(self, run_number=-1, ev_id=-1, offset = 0):
		self.particles = None
		if self.file_io is None:
			print('[e] unable to load the data because no file io is set')
			return None
		_events_match = [e for e in self.file_io.df_events if e.ev_id==ev_id and e.run_number==run_number]
		if len(_events_match) == 1:
			self.event = _events_match[0]
		else:
			print('[w] requested ev_id:', ev_id, "run_number:", run_number, 'number of matches', len(_events_match))
			return None 
		# print ('reset indexes')
		# _tmp = [p.set_user_index(0) for ip,p in enumerate(event.particles)]
		# print('loaded event:', self.current_event_in_file)
		# self.current_event_in_file = self.current_event_in_file + 1
		self.particles = fj.vectorPJ()
		for ip,p in enumerate(self.event.particles):
			p.set_user_index(offset + ip)
			self.particles.push_back(p)
		return self.particles

#random order of files; files can repeat
class DataBackgroundIO(DataIO):
	def __init__(self, **kwargs):
		self.configure_from_args(file_list='PbPb_file_list.txt', tree_name='tree_Particle')	
		super(DataBackgroundIO, self).__init__(**kwargs)

	def open_file(self):
		self.file_io = None
		self.current_event_in_file = 0
		afile = random.choice(self.list_of_files)
		print('[i] opening data file', afile)
		self.file_io = DataFileIO(file_input=afile, tree_name=self.tree_name)
		print('    number of events', self.current_file_number_of_events(), 'in tree', self.tree_name)
