#!/usr/bin/env python
import aleph
import fastjet as fj
import fjext
import fjcontrib
from tqdm import tqdm
import sys
import numpy as nd
import pandas as pd
import time
import aleph_utils
import joblib

def old_loop(e):
    parts = []
    aleph_parts = e.get_particles()
    # print(len(aleph_parts))
    for p in e.get_particles():
        psj = fj.PseudoJet(p.px(), p.py(), p.pz(), p.e())
        parts.append(psj)
    jets = jet_selector(jet_def(parts))
    #all_jets.extend(jets)

def new_loop(e):
    vparts = e.get_particles_vdoubles()
    a = nd.array(vparts)
    df = pd.DataFrame(vparts, columns=aleph.Particle.descr())
    # print()
    # print(len(df['e'].values ),   df['e'].values)
    # print(len(df['pz'].values),   df['pz'].values)
    # print(len(df['py'].values),  df['py'].values) 
    # print(len(df['px'].values),   df['px'].values)
    parts = fjext.vectorize_px_py_pz_e(df['px'].values, df['py'].values, df['pz'].values, df['e'].values)
    jets = jet_selector(jet_def(parts))
    #all_jets.extend(jets)
    # print()
    # print(df)

aleph_file="/Volumes/two/data/aleph/LEP1Data1992_recons_aftercut-001.aleph"
if len(sys.argv) > 1:
    aleph_file = sys.argv[1]

# aleph.dump(aleph_file, 2, False);
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
nev = aleph_utils.get_n_events(aleph_file)
pbar = tqdm(total=nev)

dt_stats = []
while reader.read_next_event():
    e = reader.get_event()
    aleph_parts = e.get_particles()

    pbar.update()

    start = time.time()
    old_loop(e)
    end = time.time()
    dt_old = end - start

    start = time.time()
    new_loop(e)
    end = time.time()
    dt_new = end - start

    dt_stats.append([dt_old, dt_new, dt_old/dt_new, dt_new - dt_old, len(aleph_parts)])

pbar.close()

print()
print("max speedup  : {}".format(min(dt_stats[2])))
print("min speedup  : {}".format(max(dt_stats[2])))
print("mean speedup : {}".format((max(dt_stats[2]) - min(dt_stats[2]))/len(dt_stats)))

df_timing = pd.DataFrame(dt_stats, columns=['old', 'new', 'old/new', 'new-old', 'nparts'])
joblib.dump(df_timing, 'aleph_timing.joblib')

# jet_def_lund = fj.JetDefinition(fj.cambridge_algorithm, 1.0)
# lund_gen = fjcontrib.LundGenerator(jet_def_lund)
# lunds = [lund_gen.result(j) for j in all_jets]
