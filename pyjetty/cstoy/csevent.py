# port of https://phab.hepforge.org/source/fastjetsvn/browse/contrib/contribs/ConstituentSubtractor/tags/1.4.4/example_event_wide.cc
# to python w/ heppy

import fastjet as fj
import fjcontrib

class CSEventSubtractor(object):
	def __init__(self):

		self.max_eta=4  # specify the maximal pseudorapidity for the input particles. It is used for the subtraction. Particles with eta>|max_eta| are removed and not used during the subtraction (they are not returned). The same parameter should be used for the GridMedianBackgroundEstimator as it is demonstrated in this example. If JetMedianBackgroundEstimator is used, then lower parameter should be used  (to avoid including particles outside this range). 
		self.max_eta_jet=3  # the maximal pseudorapidity for selected jets. Not used for the subtraction - just for the final output jets in this example.

		# background estimator
		self.bge_rho = fj.GridMedianBackgroundEstimator(self.max_eta, 0.2)  # Maximal pseudo-rapidity cut max_eta is used inside ConstituentSubtraction, but in GridMedianBackgroundEstimator, the range is specified by maximal rapidity cut. Therefore, it is important to apply the same pseudo-rapidity cut also for particles used for background estimation (specified by function "set_particles") and also derive the rho dependence on rapidity using this max pseudo-rapidity cut to get the correct rescaling function!  

		self.subtractor = fjcontrib.ConstituentSubtractor()  # no need to provide background estimator in this case
		self.subtractor.set_distance_type(fjcontrib.ConstituentSubtractor.deltaR)  # free parameter for the type of distance between particle i and ghost k. There  are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta  
		self.subtractor.set_max_distance(0.3)  # free parameter for the maximal allowed distance between particle i and ghost k
		self.subtractor.set_alpha(1)  # free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
		self.subtractor.set_ghost_area(0.01)  # free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
		# self.subtractor.set_do_mass_subtraction()  # use this line if also the mass term sqrt(pT^2+m^2)-pT should be corrected or not. It is necessary to specify it like this because the function set_common_bge_for_rho_and_rhom cannot be used in this case.
		self.CBS=1.0  # choose the scale for scaling the background charged particles
		self.CSS=1.0  # choose the scale for scaling the signal charged particles
		self.subtractor.set_remove_particles_with_zero_pt_and_mass(True)  # set to false if you want to have also the zero pt and mtMinuspt particles in the output. Set to true, if not. The choice has no effect on the performance. By deafult, these particles are removed - this is the recommended way since then the output contains much less particles, and therefore the next step (e.g. clustering) is faster. In this example, it is set to false to make sure that this test is successful on all systems (mac, linux).
		self.subtractor.set_grid_size_background_estimator(0.6)  # set the grid size (not area) for the background estimation with GridMedianBackgroundEstimation which is used within CS correction using charged info 

		self.subtractor.set_max_eta(self.max_eta)  # parameter for the maximal eta cut
		self.subtractor.set_background_estimator(self.bge_rho)  # specify the background estimator to estimate rho.

		print(self.subtractor.description())

	def process_event(self, parts):
		pass