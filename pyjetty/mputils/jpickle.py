from pyjetty.mputils.generic_object import GenericObject
import pickle
import sys

class VectorData(object):
	def __init__(self, pt, eta, phi, m, y, uidx=-1):
		self.set(pt, eta, phi, m, y, uidx)

	def set(self, pt, eta, phi, m, y, uidx):
		self.pt = pt
		self.eta = eta
		self.phi = phi
		self.m = m
		self.y = y
		self.uidx = uidx

	@staticmethod
	def convert_psj(j):
		return VectorData(j.pt(), j.eta(), j.phi(), j.m(), j.rapidity(), j.user_index())


class JetPicklePSJ(object):
	def __init__(self, j, constituents):
		self.j = j
		self.constituents = [p for p in constituents]

	@staticmethod
	def convert_psj(j):
		jv = VectorData(j.pt(), j.eta(), j.phi(), j.m(), j.rapidity(), j.user_index())
		cva = []
		for c in j.constituents():
			cv = VectorData(c.pt(), c.eta(), c.phi(), c.m(), c.rapidity(), c.user_index())
			cva.append(cv)
		return JetPicklePSJ(jv, cva)

class JetPickleIO(GenericObject):
	def __init__(self, **kwargs) -> None:
		super(JetPickleIO, self).__init__(**kwargs)
  
	def add_jet(self, j):
		if self.jets is None:
			self.jets = []
		self.jets.append(JetPicklePSJ.convert_psj(j))

	def write_to_file(self, fname=None):
		if fname:
			self.fname = fname
		if self.fname is None:
			self.fname = 'JetPickleIO.pkl'
		with open(self.fname, 'wb') as fout:
			pickle.dump(self.jets, fout)
			print('[i] written', self.fname, file=sys.stderr)
   
	def load_from_file(self, fname=None):
		if fname:
			self.fname = fname
		if self.fname is None:
			self.fname = 'JetPickleIO.pkl'
		with open(self.fname, 'rb') as fin:
			self.jets = pickle.load(fin)
