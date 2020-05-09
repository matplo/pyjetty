import numpy as np
import array
import fnmatch
import os
import sys

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


class ColorS(object):
	def str(*args):
		_s = ' '.join([str(s) for s in args])
		return _s
	def red(*s): 
		return '\033[91m{}\033[00m'.format(ColorS.str(*s))
	def green(*s): 
		return '\033[92m{}\033[00m'.format(ColorS.str(*s))
	def yellow(*s): 
		return '\033[93m{}\033[00m'.format(ColorS.str(*s))
	def blue(*s): 
		return '\033[34m{}\033[00m'.format(ColorS.str(*s))
	def light_purple(*s): 
		return '\033[94m{}\033[00m'.format(ColorS.str(*s))
	def purple(*s): 
		return '\033[95m{}\033[00m'.format(ColorS.str(*s))
	def cyan(*s): 
		return '\033[96m{}\033[00m'.format(ColorS.str(*s))
	def light_gray(*s): 
		return '\033[97m{}\033[00m'.format(ColorS.str(*s))
	def no_color(*s): 
		return '\033[00m{}\033[00m'.format(ColorS.str(*s))
	def black(*s):
		return '\033[98m{}\033[00m'.format(ColorS.str(*s))
	def __init__(self):
		pass
	# credit: https://www.geeksforgeeks.org/print-colors-python-terminal/
	# Python program to print 
	# colored text and background 
	def print_format_table():
		""" 
		prints table of formatted text format options 
		"""
		for style in range(8): 
			for fg in range(30, 38): 
				s1 = '' 
				for bg in range(40, 48): 
					format = ';'.join([str(style), str(fg), str(bg)]) 
					s1 += '\x1b[%sm %s \x1b[0m' % (format, format) 
				print(s1) 
			print('\n') 

def pwarning(*args, file=sys.stderr):
	print(ColorS.yellow('[w]', *args), file=file)

def pdebug(*args, file=sys.stderr):
	print(ColorS.purple('[d]', *args), file=file)

def perror(*args, file=sys.stderr):
	print(ColorS.red('[e]', *args), file=file)

def pinfo(*args, file=sys.stdout):
	print(ColorS.green('[i]', *args), file=file)

def pindent(*args, file=sys.stdout):
	print(ColorS.no_color('   ', *args), file=file)

# think about thread safe implementation
# use unique file names... for example?
class UniqueString(object):
	locked_strings = []
	def __init__(self, base=None):
		self.base = base

	def _unique(base=None):
		i = 0
		retstring = base
		if retstring is None:
			retstring = 'UniqueString_0'
		else:
			retstring = '{}_{}'.format(str(base), i)
		while retstring in UniqueString.locked_strings:
			retstring = '{}_{}'.format(retstring.split('_')[0], i)
			i = i + 1
		UniqueString.locked_strings.append(retstring)
		return retstring

	def str(self, base=None):
		if base:
			self.base = base
		return UniqueString._unique(self.base)

	def str(base=None):
		return UniqueString._unique(base)


class NoneSetWrapper(object):
	def __init__(self, name):
		self.name = name
	def __getattr__(self, key):
		try:
			return self.__dict__[key]
		except:
			print(ColorS.red('[w] {} : {} attribute is not known'.format(self.name, key)), file=sys.stderr)
			self.__setattr__(key, None)
		return None
	def description(self):
		return 'NoneSetWrapper named {}'.format(self.name)

class NoneSetWrappers(object):
	_instance = None
	def __init__(self):
		self.wrappers = {}
	def get(self, name):
		try:
			return self.wrappers[name]
		except KeyError:
			self.wrappers[name] = NoneSetWrapper(name)
		return self.wrappers[name]
	def instance():
		if NoneSetWrappers._instance is None:
			NoneSetWrappers._instance = NoneSetWrappers()
		return NoneSetWrappers._instance


def is_iterable(o):
	result = False
	try:
		tmp_iterator = iter(o)
		result = True
	except TypeError as te:
		result = False
	return result


class MPBase(object):
	_indent = 0
	def __init__(self, **kwargs):
		self.configure_from_args(name=None, args=NoneSetWrappers.instance().get(self.__class__))
		for key, value in kwargs.items():
			self.__setattr__(key, value)
		# if getattr(self, 'args'):
		# 	self.copy_attributes(self.args)
		if self.name is None:
			self.name = UniqueString.str(type(self))

	def configure_from_args(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def copy_attributes(self, ns):
		for key in ns.__dict__:
			self.__setattr__(key, getattr(ns, key))

	def __str__(self):
		s = []
		variables = self.__dict__.keys()
		MPBase._indent = MPBase._indent + 1
		_sindent = ''.join(['   ' for i in range(MPBase._indent)])
		for v in variables:
			_s = '{} ({}) = {}'.format(v, type(self.__dict__[v]), self.__dict__[v])
			_descr_method = getattr(self.__dict__[v], 'description', None)
			if callable(_descr_method):
				_s = '{} (+description) = {} {}'.format(v, type(self.__dict__[v]), _descr_method())
			if len(_s) > 500:
				if is_iterable(self.__dict__[v]):
					_s = '{} ({}) = {}'.format(v, len(self.__dict__[v]), self.__dict__[v])
			if '\n' not in _s:
				_st = (_s[:500] + '..') if len(_s) > 500 else _s
			else:
				_st = _s
			s.append(_st)
		if MPBase._indent > 1:
			sret = "..\n{} {} ({}) with \n{} -  {}".format(_sindent, self.name, '', _sindent, '\n{} -  '.format(_sindent).join(s))
		else:
			sret = "\n{}[i] {} ({}) with \n{} -  {}".format(_sindent, self.name, type(self), _sindent, '\n{} -  '.format(_sindent).join(s))
		MPBase._indent = MPBase._indent - 1
		return sret

	def description(self):
		return self.__str__()

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
