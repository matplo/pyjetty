#!/usr/bin/env python

import uproot

fname = '/Users/ploskon/data/HFtree_trains/13-06-2019/488_20190613-0256/unmerged/child_1/0001/AnalysisResults.root'
print('[i] reading from', fname)
file = uproot.open(fname)
print(file.keys)
all_ttrees = dict(file.allitems(filterclass=lambda cls: issubclass(cls, uproot.tree.TTreeMethods)))
print(all_ttrees)
tracks = all_ttrees[b'PWGHF_TreeCreator/tree_Particle;1']
print('track keys:',tracks.keys())
tracks['run_number'].show()

print('dumping some tracks')
pds_trks = tracks.pandas.df() # entrystop=10)
print(pds_trks)

print('dumping some events')
events = all_ttrees[b'PWGHF_TreeCreator/tree_event_char;1']
pds_evs = events.pandas.df()
print('number of events', len(pds_evs))
print(pds_evs[0:5])

ntrk = 0
for i, e in pds_evs.head().iterrows():
	_ntrk = int(e['n_tracks'])
	print(i, 'n_tracks =', _ntrk)
	_trks = tracks.pandas.df(entrystart=ntrk, entrystop=ntrk + _ntrk)
	for it, t in _trks.head().iterrows():
		print(it, t)
		break
	ntrk += _ntrk

for i, e in pds_evs.head().iterrows():
	iev_id = int(e['ev_id'])
	print(pds_trks.loc[pds_trks['ev_id'] == iev_id])
