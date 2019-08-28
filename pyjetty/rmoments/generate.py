#!/usr/bin/env ipython3

import fastjet as fj
import fjcontrib
import pythia8
import pythiaext
import pythiafjext
import math
from tqdm import tqdm
import joblib
from fjutils import fjpsj
from pythiautils import configuration as pyconf

def print_jets(events):
	for i, rjets in enumerate(events):
		if len(rjets):
			print(i)
			for j in rjets:
				print(j)


sconfig_pythia = ["Beams:eCM = 8000.", "HardQCD:all = on", "PhaseSpace:pTHatMin = 100."]
pythia = pyconf.create_and_init_pythia(sconfig_pythia)

# print the banner first
fj.ClusterSequence.print_banner()
print()
# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(1)
sd = fjcontrib.SoftDrop(0, 0.1, 1.0)

all_jets = []
all_jets_py = []
all_parts = []
for iEvent in tqdm(range(1000), 'event'):
	if not pythia.next():
		continue
	parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
	all_parts.append(fjpsj.pyfj_list(parts))
	jets = fj.sorted_by_pt(jet_selector(jet_def(parts)))
	ev_jets_py = []
	for idx, j in enumerate(jets):
		j.set_user_index(idx)
		jpy = fjpsj.pyfj_from_psj(j)
		ev_jets_py.append(jpy)
	all_jets_py.append(ev_jets_py)
	all_jets.extend(jets)

print('found jets:', len(all_jets), len(all_jets_py))
joblib.dump(all_jets_py, 'testing.joblib')
# print_jets(all_jets_py)
read_jets = joblib.load('testing.joblib')
# print_jets(read_jets)

joblib.dump(all_parts, 'testing_parts.joblib')
read_parts = joblib.load('testing_parts.joblib')

all_jets_py = []
for e in tqdm(read_parts, 'read event'):
	fjpjevent = []
	for p in e:
		fjpj = fj.PseudoJet(p[0], p[1], p[2], p[3])
		fjpjevent.append(fjpj)
	jets = fj.sorted_by_pt(jet_selector(jet_def(fjpjevent)))
	ev_jets_py = []
	for idx, j in enumerate(jets):
		j.set_user_index(idx)
		jpy = fjpsj.pyfj_from_psj(j)
		ev_jets_py.append(jpy)
	all_jets_py.append(ev_jets_py)
# print_jets(all_jets_py)
