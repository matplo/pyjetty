import ROOT
# https://root-forum.cern.ch/t/temporarily-turn-off-warnings-in-pyroot/18837

class Quiet(object):
	"""Context manager for silencing certain ROOT operations.  Usage:
	with Quiet(level = ROOT.kInfo+1):
	   foo_that_makes_output

	You can set a higher or lower warning level to ignore different
	kinds of messages.  After the end of indentation, the level is set
	back to what it was previously.
	"""
	def __init__(self, level=ROOT.kInfo + 1):
		self.level = level

	def __enter__(self):
		self.oldlevel = ROOT.gErrorIgnoreLevel
		ROOT.gErrorIgnoreLevel = self.level

	def __exit__(self, type, value, traceback):
		ROOT.gErrorIgnoreLevel = self.oldlevel


class QuietInfo(Quiet):
	def __init__(self):
		super().__init__(ROOT.kInfo + 1)


class QuietWarning(Quiet):
	def __init__(self):
		super().__init__(ROOT.kWarning + 1)


class QuietError(Quiet):
	def __init__(self):
		super().__init__(ROOT.kError + 1)
