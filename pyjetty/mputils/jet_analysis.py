from pyjetty.mputils.mputils import MPBase
import fastjet as fj
import fjcontrib
import fjext
import fjtools

class JetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, jet_pt_min=0.0, particles=None, explicit_ghosts=True)
		super(JetAnalysis, self).__init__(**kwargs)

		self.particle_selector = fj.SelectorAbsEtaMax(self.particle_eta_max)

		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		# self.jet_eta_max = self.particle_eta_max - 0.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(self.jet_pt_min) & fj.SelectorAbsEtaMax(self.jet_eta_max)
		if self.explicit_ghosts:
			self.jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		else:
			self.jet_area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.particle_eta_max))

		if self.particles:
			self.analyze_event(self.particles)

	def analyze_event(self, parts):
		self.particles = self.particle_selector(parts)
		if len(self.particles) < 1:
			self.cs = None
			self.jets = []
		else:
			self.cs = fj.ClusterSequenceArea(self.particles, self.jet_def, self.jet_area_def)
			self.jets = fj.sorted_by_pt(self.jet_selector(self.cs.inclusive_jets()))

	def jets_as_psj_vector(self):
		self.psj_jet_vector = fj.vectorPJ()
		# _tmp = [self.psj_jet_vector.push_back(j) for j in self.jets if not j.is_pure_ghost()]
		_tmp = [self.psj_jet_vector.push_back(j) for j in self.jets]
		return self.psj_jet_vector

class JetAnalysisPerJet(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, input_jets=None)
		super(JetAnalysisPerJet, self).__init__(**kwargs)
		self.jets = []
		self.jas = []
		for _j in self.input_jets:
			if len(_j.constituents()):
				_ja = JetAnalysis(jet_R=self.jet_R, jet_algorithm=self.jet_algorithm, particle_eta_max=self.particle_eta_max, particles=_j.constituents())
				self.jets.extend(_ja.jets)
				self.jas.append(_ja)


class JetAnalysisWithRho(JetAnalysis):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, particles=None, explicit_ghosts=True)

		self.bg_rho_range = fj.SelectorAbsEtaMax(self.particle_eta_max)
		self.bg_jet_def = fj.JetDefinition(fj.kt_algorithm, self.jet_R)
		if self.explicit_ghosts:
			self.bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		else:
			self.bg_area_def = fj.AreaDefinition(fj.active_area, fj.GhostedAreaSpec(self.particle_eta_max))
		self.bg_estimator = fj.JetMedianBackgroundEstimator(self.bg_rho_range, self.bg_jet_def, self.bg_area_def)
		self.rho = 0

		super(JetAnalysisWithRho, self).__init__(**kwargs)

	def analyze_event(self, parts):
		self.particles = self.particle_selector(parts)
		if len(self.particles) < 1:
			self.rho = 0.0
			self.cs = None
			self.jets = []
			self.corr_jet_pt = []
		else:
			self.bg_estimator.set_particles(self.particles)
			self.rho = self.bg_estimator.rho()
			self.sigma = self.bg_estimator.sigma()
			self.cs = fj.ClusterSequenceArea(self.particles, self.jet_def, self.jet_area_def)
			self.jets = fj.sorted_by_pt(self.jet_selector(self.cs.inclusive_jets()))
			self.corr_jet_pt = [j.pt() - j.area() * self.rho for j in self.jets]

	def corrected_pt(self, jet):
		return jet.pt() - jet.area() * self.rho

	def corrected_pt_plus_sigma(self, jet):
		return jet.pt() - jet.area() * (self.rho + self.sigma)

def matched_pt_constituent(c0, j1):
	_pt = [c.pt() for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index()]
	if len(_pt) == 1:
		return _pt[0]
	if len(_pt) > 1:
		print('[e] more than one match of constituents?')
		return _pt[0]
	return 0


def g_matched_pt_constituent(c0, j1):
	g = (c.pt() for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index())
	return next(g)


def matched_pt(j0, j1):
	# pts = [g_matched_pt_constituent(c0, j1) for c0 in j0.constituents()]
	pts = [c.pt() 
			for c0 in j0.constituents()
			for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index()]
	idxs = [c.user_index()
			for c0 in j0.constituents()
			for i, c in enumerate(j1.constituents()) if c.user_index() == c0.user_index()]
	pt_sum = sum(pts)
	return pt_sum / j1.pt()


def fill_tree_data(jet, tw, sd, rho, iev=None, weight=None, sigma=None):
	tw.clear()
	if iev:
		tw.fill_branch('ev_id', iev)
	if weight:
		tw.fill_branch('weight', weight)
	if sigma:
		tw.fill_branch('sigma', sigma)

	good_jet = 0.0
	if len(jet.constituents()) > 0:
		if fj.sorted_by_pt(jet.constituents())[0].pt() < 100.:
			good_jet = 1.0
	tw.fill_branch('good', good_jet)

	sd_jet = sd.result(jet)
	sd_info_jet = fjcontrib.get_SD_jet_info(sd_jet)

	tw.fill_branch('rho', rho)

	tw.fill_branch('j', jet)
	tw.fill_branch('j_ptc', jet.pt() - jet.area() * rho)
	tw.fill_branch('sd_j', sd_jet)
	tw.fill_branch('sd_j_cpt', sd_jet.pt() - sd_jet.area() * rho)
	tw.fill_branch('sd_j_z', sd_info_jet.z)
	tw.fill_branch('sd_j_dR', sd_info_jet.dR)

	pe1 = fj.PseudoJet()
	pe2 = fj.PseudoJet()
	has_parents = sd_jet.has_parents(pe1, pe2)
	tw.fill_branch('j_p1', pe1)
	tw.fill_branch('j_p2', pe2)
	if has_parents:
		tw.fill_branch('j_p1_ptc', pe1.pt() - pe1.area() * rho)
		tw.fill_branch('j_p2_ptc', pe2.pt() - pe2.area() * rho)
	else:
		tw.fill_branch('j_p1_ptc', -1000)
		tw.fill_branch('j_p2_ptc', -1000)

	tw.fill_tree()
	return jet


def fill_tree_matched(signal_jet, emb_jet, tw, sd, rho, iev=None, weight=None, sigma=None):
	tw.clear()
	mpt = fjtools.matched_pt(emb_jet, signal_jet)
	if mpt <= 0.5:
		return None

	tw.fill_branch('j_mpt', mpt)

	if iev:
		tw.fill_branch('ev_id', iev)
	if weight:
		tw.fill_branch('weight', weight)
	if sigma:
		tw.fill_branch('sigma', sigma)

	sd_signal_jet = sd.result(signal_jet)
	sd_info_signal_jet = fjcontrib.get_SD_jet_info(sd_signal_jet)
	sd_emb_jet = sd.result(emb_jet)
	sd_info_emb_jet = fjcontrib.get_SD_jet_info(sd_emb_jet)

	tw.fill_branch('rho', rho)

	tw.fill_branch('j', signal_jet)
	tw.fill_branch('sd_j', sd_signal_jet)
	tw.fill_branch('sd_j_z', sd_info_signal_jet.z)
	tw.fill_branch('sd_j_dR', sd_info_signal_jet.dR)
	tw.fill_branch('j_nc', len(signal_jet.constituents()))

	tw.fill_branch('ej', emb_jet)
	tw.fill_branch('ej_ptc', emb_jet.pt() - emb_jet.area() * rho)
	tw.fill_branch('sd_ej', sd_emb_jet)
	tw.fill_branch('sd_ej_cpt', sd_emb_jet.pt() - sd_emb_jet.area() * rho)
	tw.fill_branch('sd_ej_z', sd_info_emb_jet.z)
	tw.fill_branch('sd_ej_dR', sd_info_emb_jet.dR)

	p1 = fj.PseudoJet()
	p2 = fj.PseudoJet()
	has_parents_signal = sd_signal_jet.has_parents(p1, p2)
	# print('signal_jet:', has_parents, len(p1.constituents()), len(p2.constituents()))
	tw.fill_branch('j_p1', p1)
	tw.fill_branch('j_p2', p2)


	pe1 = fj.PseudoJet()
	pe2 = fj.PseudoJet()
	has_parents_emb = sd_emb_jet.has_parents(pe1, pe2)
	tw.fill_branch('ej_p1', pe1)
	tw.fill_branch('ej_p2', pe2)
	if has_parents_emb:
		tw.fill_branch('ej_p1_ptc', pe1.pt() - pe1.area() * rho)
		tw.fill_branch('ej_p2_ptc', pe2.pt() - pe2.area() * rho)
	else:
		tw.fill_branch('ej_p1_ptc', -1000)
		tw.fill_branch('ej_p2_ptc', -1000)

	mpt1 = -1.0 # not passed SD
	mpt2 = -1.0 # not passed SD

	if has_parents_signal and has_parents_emb:
		mpt1 = fjtools.matched_pt(pe1, p1)
		mpt2 = fjtools.matched_pt(pe2, p2)
	tw.fill_branch('mpt1', mpt1)
	tw.fill_branch('mpt2', mpt2)

		# print('signal_jet:', has_parents, len(pe1.constituents()), len(pe2.constituents()))
		# print('emb_jets', has_parents, len(pe1.constituents()), len(pe2.constituents()))

	# for c in pe2.constituents():
	# 	cp1 = fj.PseudoJet()
	# 	cp2 = fj.PseudoJet()
	# 	print(' - ', c.has_parents(cp1, cp2))

	#tw.fill_branch('jsd', sd_j)
	#tw.fill_branch('jm', ej)
	tw.fill_tree()
	return emb_jet


def sum_pt(parts):
	pts = [p.pt() for p in parts]
	return sum(pts)


def mean_pt(parts):
	return sum_pt(parts) / len(parts)


def remove_jets(parts, jets):
	psjv = fj.vectorPJ()
	reject_indexes = [p.user_index() for j in jets for p in j.constituents() if p.user_index() >= 0]
	# print ('indexes leading jets:', len(reject_indexes), sorted(reject_indexes))
	_tmp = [psjv.push_back(p) for p in parts if p.user_index() not in reject_indexes]
	return psjv
