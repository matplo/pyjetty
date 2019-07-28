#!/usr/bin/env python
import aleph
import fastjet as fj
import fjext
import fjcontrib
from tqdm import tqdm
import sys
import numpy as nd
import pandas as pd

aleph_file="/Volumes/two/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
if len(sys.argv) > 1:
    aleph_file = sys.argv[1]

aleph.dump(aleph_file, 2, False);
# aleph.dump(aleph_file, -1, True);

# print the banner first
fj.ClusterSequence.print_banner()
print()
# set up our jet definition and a jet selector
jet_R0 = 0.4
jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
#jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(200.0) & fj.SelectorAbsEtaMax(1)
jet_selector = fj.SelectorPtMin(0.0) & fj.SelectorPtMax(200.0)

all_jets = []
reader = aleph.Reader(aleph_file)
pbar = tqdm()
while reader.read_next_event():
    pbar.update()
    e = reader.get_event()
    parts = []
    # aleph_parts = e.get_particles()
    # print(len(aleph_parts))
    # for p in e.get_particles():
    #     psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
    #     parts.append(psj)
    # jets = jet_selector(jet_def(parts))
    # all_jets.extend(jets)
    vparts = e.get_particles_vdoubles()
    a = nd.array(vparts)
    df = pd.DataFrame(vparts, columns=aleph.Particle.descr())
    print()
    print(df)

# jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
# lund_gen = fjcontrib.LundGenerator(jet_def_lund)
# lunds = [lund_gen.result(j) for j in all_jets]
