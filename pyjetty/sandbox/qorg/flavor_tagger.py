from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os

import pythia8
import pythiaext
import pythiafjext


def FT_psj_to_str(p):
	psj_s = 'uidx={:3}\tpt={:7.3f}\tphi={:7.3f}\teta={:7.3f}\ty={:7.3f}'.format(p.user_index(), p.perp(), p.phi(), p.eta(), p.rap())
	py_s = 'no pythia info'
	pyp = pythiafjext.getPythia8Particle(p)
	if pyp:
		py_s = 'pt={:7.3f}\tphi={:7.3f}\teta={:7.3f}\ty={:7.3f}\t{}'.format(pyp.pT(), pyp.phi(), pyp.eta(), pyp.y(), pyp.name())
	return '\t|\t'.join([psj_s, py_s])

def FT_print_psj(jet, s='', show_constituents=False):
	if len(s) < 1:
		s = '*-*'
	print(' {} \t {}'.format(s, FT_psj_to_str(jet)))
	_stmp = ''.join([' ' for i in range(len(s))])
	if show_constituents is True:
		_ = [print('{} \t {}'.format(_stmp, FT_psj_to_str(_c))) for _c in fj.sorted_by_pt(jet.constituents())]


class FlavorTaggerUtil(object):
	def __init__(self, pythia, background_parts=None, bg_index=100000):
		self.pythia = pythia
		self.flavor_indexes = bg_index + int(1e5)
		self.background_parts = background_parts
		self.bg_index = bg_index
		self.prepare_pythia_particles()

	def _make_psj_from_iterable(self, l):
		_psj = fj.vectorPJ()
		_ = [_psj.push_back(p) for p in l]
		return _psj

	def filter_status(self, parts):
		# status 21-29: Particles from the hardest subprocess
		# status 31-39: Particles from subsequent subprocesses in multiple interactions
		_out = fj.vectorPJ()
		for p in parts:
			# note: we assume we do not care for the actual negative status
			_status = abs(pythiafjext.getPythia8Particle(p).status())
			if p.perp() <= 0:
				continue
			if (_status >= 21):  # Particles from the hardest subprocess start at -21
				_out.push_back(p)
		return _out

	def prepare_pythia_particles(self):
		_partons = pythiafjext.vectorize_select(self.pythia, [pythiafjext.kParton], 0, True)
		self.partons = self.filter_status(_partons)
		self.partons_ghosts = self._make_psj_from_iterable([_p * (1./_p.perp() * 1e-5) for _p in self.partons if _p.perp() > 0])
		_ = [self.partons_ghosts[i].set_user_index(self.partons_ghosts[i].user_index() + self.flavor_indexes) for i in range(len(self.partons_ghosts))]
		return

	# use the output of this to analyze/run jet finder ...
	def merged_particle_vector(self, signal_parts, pTcut=0.0):
		if type(signal_parts) != type(fj.vectorPJ):
			signal_parts = self._make_psj_from_iterable(signal_parts)
		_ = [signal_parts.push_back(p) for p in self.partons_ghosts]
		if self.background_parts:
			_ = [signal_parts.push_back(p) for p in self.background_parts]
		return signal_parts

	# note tagging goes with the highest pT parton within the jet as caught by the reco using partons as ghosts...
	def flavor_tag_psj(self, jet):
		_pindexes = [p.user_index() - self.flavor_indexes for p in jet.constituents() if p.user_index() >= self.flavor_indexes]
		_partons_in_the_jet = sorted([p for p in self.partons if p.user_index() in _pindexes], key=lambda x: x.user_index())
		_partons_in_the_jet = fj.sorted_by_pt(
			[p for p in self.partons if p.user_index() in _pindexes])
		if len(_partons_in_the_jet) < 1:
			return None
		lead_parton = _partons_in_the_jet[0]
		return lead_parton

	def flavor_tag_pythia_particle(self, jet):
		lead_parton = self.flavor_tag_psj(jet)
		if lead_parton:
			py_lead_parton = pythiafjext.getPythia8Particle(lead_parton)
			return py_lead_parton
		return lead_parton

	def flavor_tag_id(self, jet):
		py_lead_parton = self.flavor_tag_pythia_particle(jet)
		if py_lead_parton:
			return py_lead_parton.id()
		return py_lead_parton

	def flavor_tag_id_save(self, jet):
		_pindexes = [p.user_index() - self.flavor_indexes for p in jet.constituents() if p.user_index() >= self.flavor_indexes]
		_partons_in_the_jet = [p for p in self.partons if p.user_index() in _pindexes]
		_partons_sorted = fj.sorted_by_pt(_partons_in_the_jet)
		if len(_partons_sorted) < 1:
			return None
		lead_parton = _partons_sorted[0]
		py_lead_parton = pythiafjext.getPythia8Particle(lead_parton)
		return py_lead_parton.id()

	def is_jet_purely_parton_ghost(self, jet):
		_pindexes = [p.user_index() - self.flavor_indexes for p in jet.constituents() if p.user_index() >= self.flavor_indexes]
		if len(_pindexes) == len(jet.constituents()):
			return True
		return False
