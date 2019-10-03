import ROOT
import os


class TreeUtils(object):
	def __init__(self, tree):
		self.fname = None
		self.fin = None
		self.tree = tree
		self.bnames = []
		try:
			self.tname = self.tree.GetName()
			self.getTbranches()
		except:
			self.getT(tree)
			self.getTbranches()

	def getT(self, tpath):
		self.fname = tpath.split(".root")[0] + ".root"
		self.fin = ROOT.TFile(self.fname)
		self.tname = os.path.basename(tpath.split(".root")[1].replace(".root", ""))
		tdir = os.path.dirname(tpath.split(".root")[1].replace(".root", "")).lstrip("/")
		# print ("fname:", fname, "tdir:", tdir, "tname:", tname)
		_tpath = self.tname
		if len(tdir) == 0:
			din = self.fin
		else:
			# print ('getting ', tdir)
			_tdir = self.fin.Get(tdir)
			din = _tdir		
			# print ('din:', din)
		self.tree = din.Get(_tpath)

	def getTbranches(self):
		self.bnames = []
		if self.tree:
			lob = self.tree.GetListOfBranches()
			for b in lob:
				self.bnames.append(b.GetName())
		else:
			# print("unable to get the tree")
			print("[error] no tree...")