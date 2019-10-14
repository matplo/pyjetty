from pyjetty.mputils import MPBase
import fastjet as fj
import fjcontrib
import fjext

class JetAnalysis(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(jet_R=0.4, jet_algorithm=fj.antikt_algorithm, particle_eta_max=0.9, particles=None)
		super(JetAnalysis, self).__init__(**kwargs)

		self.bg_rho_range = fj.SelectorAbsEtaMax(self.particle_eta_max)
		self.bg_jet_def = fj.JetDefinition(fj.kt_algorithm, self.jet_R)
		self.bg_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))
		self.bg_estimator = fj.JetMedianBackgroundEstimator(self.bg_rho_range, self.bg_jet_def, self.bg_area_def)
		self.rho = 0

		self.jet_eta_max = self.particle_eta_max - self.jet_R * 1.05
		self.jet_def = fj.JetDefinition(self.jet_algorithm, self.jet_R)
		self.jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorAbsEtaMax(self.jet_eta_max)
		self.jet_area_def = fj.AreaDefinition(fj.active_area_explicit_ghosts, fj.GhostedAreaSpec(self.particle_eta_max))

		if self.particles:
			self.analyze_event(self.particles)

	def analyze_event(self, parts):
		if len(parts) < 1:
			self.rho = 0.0
			self.cs = None
			self.jets = []
			self.corr_jet_pt = []
		else:
			self.bg_estimator.set_particles(parts)
			self.rho = self.bg_estimator.rho()
			self.cs = fj.ClusterSequenceArea(parts, self.jet_def, self.jet_area_def)
			self.jets = fj.sorted_by_pt(self.jet_selector(self.cs.inclusive_jets()))
			self.corr_jet_pt = [j.pt() - j.area() * self.rho for j in self.jets]


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


def fill_tree_data(jet, tw, sd, rho, iev=None, weight=None):
	if iev:
		tw.fill_branch('ev_id', iev)
	if weight:
		tw.fill_branch('weight', weight)

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


def fill_tree_matched(signal_jet, emb_jet, tw, sd, rho, iev=None, weight=None):
	if matched_pt(emb_jet, signal_jet) <= 0.5:
		return None

	if iev:
		tw.fill_branch('ev_id', iev)
	if weight:
		tw.fill_branch('weight', weight)

	sd_signal_jet = sd.result(signal_jet)
	sd_info_signal_jet = fjcontrib.get_SD_jet_info(sd_signal_jet)
	sd_emb_jet = sd.result(emb_jet)
	sd_info_emb_jet = fjcontrib.get_SD_jet_info(sd_emb_jet)

	tw.fill_branch('rho', rho)

	tw.fill_branch('j', signal_jet)
	tw.fill_branch('sd_j', sd_signal_jet)
	tw.fill_branch('sd_j_z', sd_info_signal_jet.z)
	tw.fill_branch('sd_j_dR', sd_info_signal_jet.dR)

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
		mpt1 = matched_pt(pe1, p1)
		mpt2 = matched_pt(pe2, p2)
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

