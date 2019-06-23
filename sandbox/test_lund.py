#!/usr/bin/env python3
import fastjet as fj
import pythia8
from recursivetools import pyrecursivetools as rt
from lundplane import pylundplane as lund
from pythiafjtools import pypythiafjtools as pyfj
import math
from tqdm import tqdm


def deltas(jets, jets0):
    for i in range(len(jets)):
        yield jets0[i].perp() - jets[i].perp()


def create_and_init_pythia(config_strings=[]):
    pythia = pythia8.Pythia()
    for s in config_strings:
        pythia.readString(s)
    for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0"]:
        pythia.readString(extra_s)
    if pythia.init():
        return pythia
    return None


sconfig_pythia = ["Beams:eCM = 8000.", "HardQCD:all = on", "PhaseSpace:pTHatMin = 100."]
pythia = create_and_init_pythia(sconfig_pythia)

# print the banner first
fj.ClusterSequence.print_banner()
print()
# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
jet_selector = fj.SelectorPtMin(100.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(1)
sd = rt.SoftDrop(0, 0.1, 1.0)

all_jets = []
for iEvent in tqdm(range(10), 'event'):
    if not pythia.next(): continue
    parts = pyfj.vectorize(pythia, True, -1, 1, False)
    jets = jet_selector(jet_def(parts))
    all_jets.extend(jets)

all_sd_jets = [sd.result(j) for j in all_jets]
pts = [j.pt() for j in all_jets]
sd_pts = [j.pt() for j in all_sd_jets]
sd_delta_pt = [delta for delta in deltas(all_jets, all_sd_jets)]
nangs0 = [pyfj.angularity(j, 0.) for j in all_jets]

jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
lund_gen = lund.LundGenerator(jet_def_lund)

lunds = [lund_gen.result(j) for j in all_jets]
print('listing lund plane points... Delta, kt - for {} selected jets'.format(len(all_jets)))
for l in lunds:
    print ('- jet pT={0:5.2f} eta={1:5.2f}'.format(l[0].pair().perp(), l[0].pair().eta()))
    print ('  Deltas={}'.format([s.Delta() for s in l]))
    print ('  kts={}'.format([s.Delta() for s in l]))
