#!/usr/bin/env python3
import fastjet as fj
import pythia8
from pythiafjtools import pypythiafjtools as pyfj
import math
from tqdm import tqdm


def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	for s in config_strings:
		pythia.readString(s)
	for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0"]:
		pythia.readString(extra_s)
	if pythia.init():
		return pythia
	return None


sconfig_pythia = [	"Beams:eCM = 5000.", "HardQCD:all = on",
					"PartonLevel:ISR = off",
					"PartonLevel:MPI = off",
					"PhaseSpace:bias2Selection=on",
					"PhaseSpace:bias2SelectionPow=4",
					"PhaseSpace:bias2SelectionRef=50"]
pythia = create_and_init_pythia(sconfig_pythia)

# print the banner first
fj.ClusterSequence.print_banner()
print()
# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(1)

all_jets = []
for iEvent in tqdm(range(1000), 'event'):
	if not pythia.next(): continue
	parts = pyfj.vectorize(pythia, True, -1, 1, False)
	jets = jet_selector(jet_def(parts))
	all_jets.extend(jets)

