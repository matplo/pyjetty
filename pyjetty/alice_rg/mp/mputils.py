import numpy as np
import array

def logbins(xmin, xmax, nbins):
	if xmin <= 0:
		xmin = 1e-2
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr


class MPBase(object):
	def __init__(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def configure_constants(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def __str__(self):
		s = []
		variables = self.__dict__.keys()
		for v in variables:
			s.append('{} = {}'.format(v, self.__dict__[v]))
		return "[i] {} with \n .  {}".format(self.__class__.__name__, '\n .  '.join(s))

import sys

class CursorSpin(object):
	_cpos = 0
	_cursor = '\\|/-'
	def __init__(self):
		sys.stdout.write(' {}\r'.format(CursorSpin._cursor[CursorSpin._cpos]))
		sys.stdout.flush()
		CursorSpin._cpos = CursorSpin._cpos + 1
		if CursorSpin._cpos >= len(CursorSpin._cursor):
			CursorSpin._cpos = 0