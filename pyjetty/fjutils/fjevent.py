import fastjet as fj
import pandas as pd

class KineSelectorFactory(object):
	def __init__(self, **kwargs):
		self.ptmin = 0.0
		self.ptmax = +1.e6
		self.etamin = -1.e6
		self.etamax = +1.e6
		self.absetamax = None
		for key, value in kwargs.items():
			if key == "ptmin": self.ptmin = value
			if key == "ptmax": self.ptmax = value
			if key == "etamin": self.etamin = value
			if key == "etamax": self.etamax = value
			if key == "absetamax": self.absetamax = value
		if self.absetamax is None:
			if abs(self.etamin) == self.etamax:
				self.absetamax = self.etamax
		if self.absetamax:
			self.selector = fj.SelectorPtRange(self.ptmin, self.ptmax) & fj.SelectorAbsEtaMax(self.absetamax)
		else:
			self.selector = fj.SelectorPtRange(self.ptmin, self.ptmax) & fj.SelectorEtaMin(self.etamin) & fj.SelectorEtaMax(self.etamax)


class FJEvent(object):
	df_columns = ['evid', 'pt', 'eta', 'phi', 'area', 'ptsub', 'name']
	def __init__(self, **kwargs):
		# eventid = -1, R=0.4, algorithm=fj.antikt_algorithm, name = None):
		self.particles = []
		self.R = 0.4
		self.algorithm = fj.antikt_algorithm
		self.csaa = None
		self.jet_def = None
		self.jet_area_def = None
		self.cs = None
		self.csaa = None
		self.particle_selector = None
		self.jet_selector = None
		self.name = None
		self.id = -1
		self.ngrid = None
		for key, value in kwargs.items():
			if key == "R": self.R = value
			if key == "algorithm": self.algorithm = value
			if key == "particle_selector": self.particle_selector = value
			if key == "jet_selector": self.jet_selector = value
			if key == "name": self.name = value
			if key == "id": self.id = value
			if key == "particles": self.particles = value
		if self.name is None:
			self.name = '{}_{}'.format(self.R, self.algorithm)
		self.background_estimator = None
		self.jets_df = None

	def run_jet_finder_csaa(self):
		self.jet_def = fj.JetDefinition(self.algorithm, self.R)
		self.jet_area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.particle_selector.absetamax))
		self.csaa = fj.ClusterSequenceArea(self.particle_selector.selector(self.particles), self.jet_def, self.jet_area_def)
		self.inclusive_jets = fj.sorted_by_pt(self.csaa.inclusive_jets())
		if self.jet_selector:
			self.inclusive_jets = self.jet_selector.selector(self.inclusive_jets)
		if self.ngrid is None:
			self.ngrid = self.particle_selector.absetamax / 0.05
		_grid_spacing = self.particle_selector.absetamax / self.ngrid
		self.background_estimator = fj.GridMedianBackgroundEstimator(self.particle_selector.absetamax, _grid_spacing)
		self.background_estimator.set_particles(self.particle_selector.selector(self.particles))
		self.jets_df = pd.DataFrame([[	self.id, j.perp(), j.eta(), j.phi(), j.area(), 
										j.perp() - self.background_estimator.rho() * j.area(), self.name] 
										for j in self.inclusive_jets], 
									columns=self.df_columns)

	def leading_pt_particle(self, ptmin):
		return fj.sorted_by_pt(self.particles)[0]

	def leading_particle_select(self, selector):
		return fj.sorted_by_pt(selector(self.particles))[0]

	def leading_pt_jet(self, ptmin):
		return fj.sorted_by_pt(self.inclusive_jets)[0]

	def leading_jet_select(self, selector):
		return fj.sorted_by_pt(selector(self.inclusive_jets))[0]

