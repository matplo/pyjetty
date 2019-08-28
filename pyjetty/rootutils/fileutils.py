import uproot

def list_trees_dict(fname):
	file = uproot.open(fname)
	all_ttrees = dict(file.allitems(filterclass=lambda cls: issubclass(cls, uproot.tree.TTreeMethods)))
	return all_ttrees
