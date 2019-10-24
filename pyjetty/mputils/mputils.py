import numpy as np
import array
import fnmatch
import os

def find_files(rootdir='.', pattern='*'):
    return [os.path.join(rootdir, filename)
            for rootdir, dirnames, filenames in os.walk(rootdir)
            for filename in filenames
            if fnmatch.fnmatch(filename, pattern)]

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

	def configure_from_args(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def __str__(self):
		s = []
		variables = self.__dict__.keys()
		for v in variables:
			s.append('{} = {}'.format(v, self.__dict__[v]))
		return "\n[i] {} with \n .  {}".format(self.__class__.__name__, '\n .  '.join(s))

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

class Type(object):
	_float = type(0.)
	_int = type(0)
	_list = type([0,1])
	_tuple = type((0,1))
	def __init__(self):
		pass
	def is_float(x):
		return (float == type(x))
	def is_int(x):
		return (int == type(x))
	def is_list(x):
		return (list == type(x))
	def is_tuple(x):
		return (tuple == type(x))
	def is_dict(x):
		return (dict == type(x))
