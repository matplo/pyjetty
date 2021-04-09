import uproot

def list_trees_dict(fname):
	if uproot.version.version_info[0] == '3':
		file = uproot.open(fname)
		all_ttrees = dict(file.allitems(filterclass=lambda cls: issubclass(cls, uproot.tree.TTreeMethods)))
		return all_ttrees
	if uproot.version.version_info[0] == '4':
		with uproot.open(fname) as f:
			all_ttrees = dict(f.classnames(filter_classname="TNtuple"))
			all_ttrees.update(dict(f.classnames(filter_classname="TTree")))
			return all_ttrees
	return None
