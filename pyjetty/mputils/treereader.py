import ROOT
from pyjetty.mputils import MPBase


class RTreeReader(MPBase):
	errors = []
	def __init__(self, **kwargs):
		self.configure_from_args(	tree=None, 
									tree_name=None,
									name="RTreeWriter", 
									file_name="RTreeWriter.root", 
									fin=None,
									quiet=True,
									branches = [])
		super(RTreeReader, self).__init__(**kwargs)
		if self.tree is None:
			if self.fin is None:
				if not self.quiet:
					print('[i] new file {}'.format(ROOT.gSystem.ExpandPathName(self.file_name)))
				self.fin = ROOT.TFile(self.file_name)
				self.fin.cd()
			if self.tree_name is None:
				self.tree_name = 't'+self.name
			self.tree = self.fin.Get(self.tree_name)
		self.branch_containers = {}
		self.bad_tree_or_branch = False
		for bname in self.branches:
			self.read_branch(bname)

	def read_branch(self, bname):
		b = self.tree.GetBranch(bname)
		if not b:
			if not self.quiet:
				print('[i] RTreeReader {} tree {}: branch [{}] not found'.format(self.name, self.tree.GetName(), bname))
				if self.fin:
					print('    current file: {}'.format(self.fin.GetName()))
			else:
				self.errors.append('[i] RTreeReader {} tree {}: branch [{}] not found'.format(self.name, self.tree.GetName(), bname))
				if self.fin:
					self.errors.append('    current file: {}'.format(self.fin.GetName()))
			self.bad_tree_or_branch = True
		else:
			self.branch_containers[bname] = ROOT.std.vector('float')()
			self.tree.SetBranchAddress(bname, ROOT.AddressOf(self.branch_containers[bname]))
			setattr(self, bname, self.branch_containers[bname])

	def next_event(self):
		if self.bad_tree_or_branch:
			return False
		for i in range(self.tree.GetEntries()):
			self.tree.GetEntry(i)
			yield True
		yield False

	def close(self):
		self.fin.Close()

	def print_errors():
		for e in RTreeReader.errors:
			print(e)

def test():
	tr = RTreeReader(	tree_name='t', 
						branches = ['j_pt', 'ej_pt'],
						file_name='$HOME/devel/pyjetty/pyjetty/cstoy/output_alpha_0_dRmax_0.4_SDzcut_0.2_emb.root')
	for i in range(tr.tree.GetEntries()):
		tr.tree.GetEntry(i)
		print (tr.j_pt.size(), tr.ej_pt.size())
		for ie in range(tr.j_pt.size()):
			print(tr.j_pt[ie])


if __name__ == '__main__':
	test()