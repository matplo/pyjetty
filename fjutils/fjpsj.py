import fastjet as fj
import copy


class PSJVector(object):

	def __init__(self, px=0, py=0, pz=0, e=0):
		self.p = [0, 0, 0, 0]
		self.p[0] = px
		self.p[1] = py
		self.p[2] = pz
		self.p[3] = e
		self.constituents = []
		self.user_index = -1

	def __getitem__(self, i):
		return self.p[i]

	def __str__(self):
		return self.__repr__()

	def __repr__(self):
		_s = [' - jet with user_index={0:4d}'.format(self.user_index)]
		_s.append('   px={0:7.2f} py={1:7.2f} pz={2:7.2f} e={3:7.2f}'.format(self.p[0], self.p[1], self.p[2], self.p[3]))
		_s.append('   n constituents={0:4d}'.format(len(self.constituents)))
		return '\n'.join(_s)


class Container(object):
	def __init__(self):
		pass

	@classmethod
	def __getattr__(self, name):
		self.__setattr__(name, 0)

	@classmethod
	def __setattr__(self, name, val):
		setattr(self, name, val)

	@classmethod
	def from_kwargs(cls, **kwargs):
		obj = cls()
		for (field, value) in kwargs.items():
			setattr(obj, field, value)
		return obj


def pyfj_from_psj(fjpsj):
	j = PSJVector(px=fjpsj.px(), py=fjpsj.py(), pz=fjpsj.pz(), e=fjpsj.e())
	j.user_index = fjpsj.user_index()
	if fjpsj.has_constituents() and len(fjpsj.constituents()) > 1:
		for c in fjpsj.constituents():
			_c = pyfj_from_psj(c)
			_c.user_index = c.user_index()
			j.constituents.append(copy.deepcopy(_c))
	return j

def pyfj_list(l):
	retv = []
	for vx in l:
		nvx = pyfj_from_psj(vx)
		retv.append(nvx)
	return retv
