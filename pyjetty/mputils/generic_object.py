class GenericObject(object):
	max_chars = 1000
	def __init__(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)
		if self.args:
			if type(self.args) == dict:
				self.configure_from_dict(self.args, ignore_none=True)
			else:
				self.configure_from_dict(self.args.__dict__)

	def configure_from_args(self, **kwargs):
		for key, value in kwargs.items():
			self.__setattr__(key, value)

	def configure_from_dict(self, d, ignore_none=False):
		for k in d:
			if d[k] is None:
				continue
			self.__setattr__(k, d[k])

	def __getattr__(self, key):
		try:
			return self.__dict__[key]
		except:
			pass
		self.__setattr__(key, None)
		return self.__getattr__(key)

	def __str__(self) -> str:
		s = []
		s.append('[i] {} ({})'.format(str(self.__class__).split('.')[1].split('\'')[0], id(self)))
		for a in self.__dict__:
			if a[0] == '_':
				continue
			sval = str(getattr(self, a))
			if len(sval) > self.max_chars:
				sval = sval[:self.max_chars-4] + '...'
			s.append('   {} = {}'.format(str(a), sval))
		return '\n'.join(s)

	def __repr__(self) -> str:
		return self.__str__()

	def __getitem__(self, key):
		return self.__getattr__(key)

	def __iter__(self):
		_props = [a for a in self.__dict__ if a[0] != '_']
		return iter(_props)
