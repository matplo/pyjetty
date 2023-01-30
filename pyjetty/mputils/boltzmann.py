import ROOT
import fastjet as fj
import numpy as np
import array
from pyjetty.mputils.mputils import MPBase, UniqueString, logbins

class BoltzmannEvent(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(mean_pt=0.7, multiplicity=1, max_eta=1, max_pt=100, min_pt=0.15)
		super(BoltzmannEvent, self).__init__(**kwargs)
		if self.min_pt < 0:
			self.min_pt = 0
		self.particles = fj.vectorPJ()
		self.funbg = ROOT.TF1("funbg", "2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", self.min_pt, self.max_pt, 1);
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
		# print (self)

	def write(self):
		self.histogram_phi.Write()
		self.histogram_pt.Write()
		self.histogram_eta.Write()

	def _boltzmann(self, pt):
		return 2. / self.mean_pt * pt * math.exp(-(2. / self.mean_pt) * pt);

	def generate(self, multiplicity=None, offset=0):
		if multiplicity:
			self.multiplicity = multiplicity
		self.particles.clear()
		for n in range(0, int(self.multiplicity)):
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


class BoltzmannSubtractor(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(max_pt_subtract=1000)
		super(BoltzmannSubtractor, self).__init__(**kwargs)
		fname = UniqueString.str("BoltzmannSubtractor_function")
		# self.funbg = ROOT.TF1(fname, "[1] * 2. / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, self.max_pt_subtract, 2)
		self.funbg = ROOT.TF1(fname, "[1] / [0] * x * TMath::Exp(-(2. / [0]) * x)", 0, self.max_pt_subtract, 2)
		self.funbg.SetParameter(0, 0.7)
		self.funbg.SetParameter(1, 1)
		# npx = int(self.max_pt_subtract/0.1)
		# if npx < 100:
		# 	npx = 100
		# self.funbg.SetNpx(npx)
		# hname = UniqueString.str("BoltzmannSubtractor_histogram")
		# self.h = ROOT.TH1F(hname, hname, 2*npx, -self.max_pt_subtract, self.max_pt_subtract)
		# self.h.SetDirectory(0)

	def do_subtraction(self, parts):
		self.sparts = []
		for p in parts:
			if p.pt() > self.max_pt_subtract:
				self.sparts.append(p)
				continue
			else:
				fraction = self.funbg.Eval(p.pt()) * self.total_pt
				if p.pt() - fraction > 0:
					newp = fj.PseudoJet()
					newp.reset_PtYPhiM(p.pt() - fraction, p.rap(), p.phi(), p.m())
					self.sparts.append(newp)

	def subtracted_particles(self, parts):
		# self.h.Reset()
		_parts_pt = [p.pt() for p in parts if p.pt() < self.max_pt_subtract]
		self.total_pt = sum(_parts_pt)
		# _tmp = [self.h.Fill(pt) for pt in _parts_pt]
		# self.h.Scale(1./len(_parts_pt))
		# self.h.Fit(self.funbg, "RMNQ0")
		# print ('par 1', self.funbg.GetParameter(1))
		# print ('par 0', self.funbg.GetParameter(0))
		mult = len(_parts_pt)
		mean_pt = self.total_pt / mult
		self.funbg.SetParameter(0, mean_pt)
		self.funbg.SetParameter(1, mean_pt * 2.)
		# print ('par 1', self.funbg.GetParameter(1))
		# print ('par 0', self.funbg.GetParameter(0))
		self.do_subtraction(parts)
		return self.sparts
