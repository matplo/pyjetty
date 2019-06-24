#!/usr/bin/env ipython

import fastjet as fj
import pythia8
from pythiafjtools import pypythiafjtools as pyfj
from lundplane import pylundplane as lund
from mptools import pymptools as mpt
from tqdm import tqdm
import math
import tqdm
import os
import re
import numpy as np
import pandas as pd
import joblib

def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	for s in config_strings:
		pythia.readString(s)
	for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0", "Next:numberCount = 0"]:
		pythia.readString(extra_s)
	if pythia.init():
		return pythia
	return None


sconfig_pythia_w = ["Beams:eCM = 5000.", "HardQCD:all = on",
					"PartonLevel:ISR = off",
					"PartonLevel:MPI = off",
					"PhaseSpace:bias2Selection=on",
					"PhaseSpace:bias2SelectionPow=4",
					"PhaseSpace:bias2SelectionRef=50"]
sconfig_pythia = ["Beams:eCM = 5000.", "HardQCD:all = on",
					"PartonLevel:ISR = off",
					"PartonLevel:MPI = off",
					"PhaseSpace:pThatMin = 20"]
pythia = create_and_init_pythia(sconfig_pythia)

# print the banner first
fj.ClusterSequence.print_banner()
print()
# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
jet_selector = fj.SelectorPtMin(10.0) & fj.SelectorAbsEtaMax(2 - jet_R0 * 1.05)

all_jets = []
x1 = []
x2 = []
id1 = []
id2 = []
pthat = []
pts = []
etas = []
ijet = 0

for iEvent in tqdm.tqdm(range(100000), 'event'):
	if not pythia.next():
		continue
	parts = pyfj.vectorize(pythia, True, -10, 10, False)
	jets = jet_selector(jet_def(parts))
	all_jets.extend(jets)
	_x1 = pythia.info.x1()
	_x2 = pythia.info.x2()
	_id1 = pythia.info.id1()
	_id2 = pythia.info.id2()
	_pthat = pythia.info.pTHat()
	x1.extend([_x1 for j in jets])
	x2.extend([_x2 for j in jets])
	id1.extend([_id1 for j in jets])
	id2.extend([_id2 for j in jets])
	pthat.extend([_pthat for j in jets])
	etas.extend([j.eta() for j in jets])
	pts.extend([j.pt() for j in jets])

jpd = pd.DataFrame({'pt': pts, 'eta': etas, 'x1': x1, 'x2': x2, 'id1': id1, 'id2': id2, 'pthat': pthat})

foutname = 'jets_x1x2.data'
joblib.dump(jpd, foutname, compress=9)

print('[i] written', foutname)
