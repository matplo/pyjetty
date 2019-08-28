import fastjet as fj
import pandas as pd


class FJEvent(object):
	def __init__(self, fjpsjs = []):
		self.particles = fjpsjs
		self.particle_selector = None
		self.jet_selector = None
		self.particle_pt_min = 0
		self.particle_pt_max = 1e6
		self.particle_abseta_max = 20.
		self.jet_pt_min = 0
		self.jet_pt_max = 1e6
		self.jet_abseta_max = 20.
		self.grid_spacing = 20
		self.background_estimator = None

		# dictionaties!
		self.jet_def = {}
		self.jet_area_def = {}
		self.csaa = {}
		self.inclusive_jets_csaa = {}
		self.Rs = []

	def __getitem__(self, i):
		return self.particles[i]

	def setup_particle_selector(self, ptmin = None, ptmax = None, absetamax = None):
		if ptmin:
			self.particle_pt_min = ptmin
		if ptmax:
			self.particle_pt_max = ptmax
		if absetamax:
			self.particle_abseta_max = absetamax
		self.particle_selector = fj.SelectorPtMin(self.particle_pt_min) & fj.SelectorPtMax(self.particle_pt_max) & fj.SelectorAbsEtaMax(self.particle_abseta_max)

	def leading_pt(self, ptmin):
		return fj.sorted_by_pt(self.particles)[0]

	def leading_select(self, ptmin, ptmax, absetamax):
		_selector = fj.SelectorPtMin(ptmin) & fj.SelectorPtMax(ptmax) & fj.SelectorAbsEtaMax(absetamax)
		return fj.sorted_by_pt(_selector(self.particles))[0]

	def setup_jet_selector(self, ptmin = None, ptmax = None, absetamax = None):
		if ptmin:
			self.jet_pt_min = ptmin
		if ptmax:
			self.jet_pt_max = ptmax
		if absetamax:
			self.jet_abseta_max = absetamax
		self.jet_selector = fj.SelectorPtMin(self.jet_pt_min) & fj.SelectorPtMax(self.jet_pt_max) & fj.SelectorAbsEtaMax(self.jet_abseta_max)

	def jets_csaa(self, R, algorithm = fj.antikt_algorithm):
		maxrap = self.particle_abseta_max
		self.jet_Rs.append(R)
		self.jet_def[R] = fj.JetDefinition(algorithm, R)
		self.jet_area_def[R] = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(maxrap))
		self.csaa[R] = fj.ClusterSequenceArea(_tpsj, self.jet_def, self.jet_area_def)
		self.inclusive_jets_csaa[R] = self.csaa.inclusive_jets()
		return self.inclusive_jets_csaa[R]

	def setup_background_estimator(self, grid_spacing = None):
		if grid_spacing = None:
			grid_spacing = self.particle_abseta_max/self.grid_spacing
		self.background_estimator = fj.GridMedianBackgroundEstimator(self.particle_abseta_max, grid_spacing)

	def run_background_estimation(self):
		if self.background_estimator is None:
			self.setup_background_estimator()
		self.background_estimator.set_particles(self.particles)
		self.rho = self.background_estimator.rho()
		self.sigma = self.background_estimator.sigma()

	def jets_csaa_DF(self, R, eventid = -1):
		output_columns = ['evid', 'pt', 'eta', 'phi', 'area', 'ptsub', 'R']
		_jets_df = pd.DataFrame([[eventid, j.perp(), j.eta(), j.phi(), j.area(), j.perp() - gmbge.rho() * j.area(), R] for j in self.inclusive_jets_csaa[R]], 
								columns=output_columns)
		return _jets_df

