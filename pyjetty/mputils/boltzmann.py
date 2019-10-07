import ROOT
import fastjet as fj
import numpy as np
import array
from pyjetty.mputils import MPBase
from pyjetty.mputils import logbins

class BoltzmannEvent(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(mean_pt=0.7, multiplicity=1, max_eta=1, max_pt=100)
		super(BoltzmannEvent, self).__init__(**kwargs)
		self.particles = fj.vectorPJ()
		self.funbg = ROOT.TF1("funbg", "2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, self.max_pt, 1);
		self.funbg.SetParameter(0, self.mean_pt)
		self.funbg.SetNpx(1000)
		self.ROOT_random = ROOT.TRandom()
		self.histogram_pt = ROOT.TH1F("BoltzmannEvent_pt", "BoltzmannEvent_pt;p_{T} (GeV/c)", 100, logbins(1e-1, self.max_pt, 100))
		self.histogram_pt.SetDirectory(0)
		self.histogram_eta = ROOT.TH1F("BoltzmannEvent_eta", "BoltzmannEvent_eta;#eta", 100, -self.max_eta, self.max_eta)
		self.histogram_eta.SetDirectory(0)
		self.histogram_phi = ROOT.TH1F("BoltzmannEvent_phi", "BoltzmannEvent_phi;#varphi (rad)", 100, -ROOT.TMath.Pi(), ROOT.TMath.Pi())
		self.histogram_phi.SetDirectory(0)
		self.nEvent = 0
		print (self)

	def _boltzmann(self, pt):
		return 2. / self.mean_pt * pt * math.exp(-(2. / self.mean_pt) * pt);

	def generate(self, multiplicity=None, offset=0):
		if multiplicity:
			self.multiplicity = multiplicity
		self.particles.clear()
		for n in range(0, self.multiplicity):
			_pt  = self.funbg.GetRandom(0, self.max_pt)
			_eta = self.ROOT_random.Rndm() * self.max_eta * 2. - self.max_eta;
			_phi = self.ROOT_random.Rndm() * ROOT.TMath.Pi() * 2. - ROOT.TMath.Pi();
			_p = fj.PseudoJet()
			_p.reset_PtYPhiM (_pt, _eta, _phi, 0.0)
			_p.set_user_index(n + offset)
			self.particles.push_back(_p)
			self.histogram_pt.Fill(_pt)
			self.histogram_eta.Fill(_eta)
			self.histogram_phi.Fill(_phi)
		self.nEvent = self.nEvent + 1
		return self.particles
