#!/usr/bin/env python3
from mptools import pymptools as mpt
import fastjet as fj
from lundplane import pylundplane as lund
from tqdm import tqdm
import sys
import os

aleph_file = "/Volumes/two/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
if len(sys.argv) > 1:
    aleph_file = sys.argv[1]
if not os.path.isfile(aleph_file):
    aleph_file = "/Volumes/mp256s/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
if not os.path.isfile(aleph_file):
    aleph_file = "/Users/ploskon/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
mpt.dump(aleph_file, 2, False)
# mpt.dump(aleph_file, -1, True);

# print the banner first
fj.ClusterSequence.print_banner()
print()
# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
#jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(1)
jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(200.0)

all_jets = []
reader = mpt.Reader(aleph_file)
pbar = tqdm()
while reader.read_next_event():
    pbar.update()
    e = reader.get_event()
    parts = []
    # aleph_parts = e.get_particles()
    # print(len(aleph_parts))
    for p in e.get_particles():
        psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
        parts.append(psj)
    jets = jet_selector(jet_def(parts))
    all_jets.extend(jets)

jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
lund_gen = lund.LundGenerator(jet_def_lund)
lunds = [lund_gen.result(j) for j in all_jets]
