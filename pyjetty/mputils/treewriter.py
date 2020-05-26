import ROOT
import fastjet as fj
import fjcontrib
from pyjetty.mputils import MPBase, pwarning


def get_LundDeclusteringType():
	j = fj.PseudoJet()
	ld = fjcontrib.LundDeclustering(j, j, j)
	return type(ld)


class RTreeWriter(MPBase):
	_fj_psj_type = type(fj.PseudoJet())
	_fj_psj_vector_type = type(fj.vectorPJ())
	_fj_LundDeclustering_type = get_LundDeclusteringType()
	# _fj_sdinfo = type(fjcontrib.SDinfo())
	def __init__(self, **kwargs):
		self.configure_from_args(	tree=None, 
									tree_name=None,
									name="RTreeWriter", 
									file_name="RTreeWriter.root", 
									fout=None)
		super(RTreeWriter, self).__init__(**kwargs)
		self._warnings = []
		if self.tree is None:
			if self.fout is None:
				print('[i] new file {}'.format(self.file_name))
				self.fout = ROOT.TFile(self.file_name, 'recreate')
				self.fout.cd()
			else:
				self.name = self.fout.GetName()
				self.file_name = self.name
				self.fout.cd()
			if self.tree_name is None:
				self.tree_name = 't'+self.name
			self.tree = ROOT.TTree(self.tree_name, self.tree_name)
		self.branch_containers = {}

	def add_warning(self, s):
		if s not in self._warnings:
			self._warnings.append(s)

	def _fill_branch(self, bname, value):
		b = self.tree.GetBranch(bname)
		if not b:
			print('[i] RTreeWriter {} tree {}: creating branch [{}]'.format(self.name, self.tree.GetName(), bname))
			self.branch_containers[bname] = ROOT.std.vector('float')()
			b = self.tree.Branch(bname, self.branch_containers[bname])
		if b:
			# print('filling branch:', bname, 'at', b)
			self.branch_containers[bname].push_back(value)

	def fill_branches_attribs(self, o, attr_list=[], prefix=''):
		if len(attr_list) == 0:
			attr_list = o.__dict__
		for a in attr_list:
			self.fill_branch(prefix+a, getattr(o, a))

	def fill_branches(self, **kwargs):
		for a in kwargs:
			self.fill_branch(bname=a, value=kwargs[a])

	def fill_branch(self, bname, value, do_enumerate=False):
		# print("FILL:", self.tree_name, bname, value)
		if float == type(value) or int == type(value):
			self._fill_branch(bname, value)
			return
		if type(value) in [tuple, list, self._fj_psj_vector_type]:
			if do_enumerate:
				r = [self.fill_branch('{}_{}'.format(bname, i), x) for i,x in enumerate(value)]
			else:
				r = [self.fill_branch(bname, x) for x in value]
			return
		if dict == type(value):
			r = [self.fill_branch('{}_{}'.format(bname, i), x) for i, x in value.items()]
			return
		if self._fj_psj_type == type(value):
			if value.has_area():
				self.fill_branch(bname, {	'pt' 	: value.pt(), 
											'phi' 	: value.phi(), 
											'eta' 	: value.eta(),
											'a' 	: value.area()})
			else:
				self.fill_branch(bname, {'pt' : value.pt(), 'phi' : value.phi(), 'eta' : value.eta()})
			return
		if self._fj_LundDeclustering_type == type(value):
			self.fill_branch(bname, { 
										'm' : value.m(),
										'z' : value.z(),
										'Delta' : value.Delta(),
										'kt' : value.kt(),
										'kappa' : value.kappa(),
										'psi' : value.psi(),
										'p' : value.pair(),
										's1' : value.harder(),
										's2' : value.softer(),
										'tf' : value.z() * value.Delta() * value.Delta()
									})
			return
		if bool == type(value):
			self._fill_branch(bname, value)
			return
		try:
			_val = float(value)
			self._fill_branch(bname, _val)
			self.add_warning('converted {} to float for branch {}'.format(type(value), bname))
			return			
		except:
			pass
		self.add_warning('do not know how to fill tree {} branch {} for type {} - ignored'.format(self.tree_name, bname, type(value)))

	def clear(self):
		for k in self.branch_containers:
			self.branch_containers[k].clear()

	def fill_tree(self):
		self.tree.Fill()
		self.clear()

	def write_and_close(self):
		print('[i] writing {}'.format(self.fout.GetName()))
		self.fout.Write()
		self.fout.Purge()
		self.fout.Close()

	def __del__(self):
		for w in self._warnings:
			pwarning(self.tree_name, ':', w)

def example():
	tw = RTreeWriter()
	print(tw)
	tw.fill_branch('b', 10)
	tw.fill_branch('b', 12.)
	tw.fill_branch('bl', [1, 2, 3], do_enumerate=True)
	tw.fill_branch('bt', (10, 20, 30.))
	psj = fj.PseudoJet()
	tw.fill_branch('jet', psj)
	tw.fill_branch('jet', psj)

	v = fj.vectorPJ()
	_v = fj.PseudoJet(1,2,3,4)
	v.push_back(_v)
	v.push_back(_v)
	v.push_back(_v)

	tw.fill_branch('jets', v)

	tw.fill_branch('bd', {'x':10, 'y':20, 'z':30.})
	tw.fill_tree()
	tw.write_and_close()


if __name__ == '__main__':
	example()