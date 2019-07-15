#!/usr/bin/env python

import ROOT as r
import pandas as pd
fname = '/Users/ploskon/data/HFtree_trains/13-06-2019/488_20190613-0256/unmerged/child_1/0001/AnalysisResults.root'
print('[i] reading from', fname)
tdf_parts = r.RDF.MakeRootDataFrame("PWGHF_TreeCreator/tree_Particle", fname)
np_parts = tdf_parts.AsNumpy()
df_parts = pd.DataFrame(np_parts)
df_parts_sample = df_parts.sample(frac=0.1)
print(df_parts_sample)

tdf_ev = r.RDF.MakeRootDataFrame("PWGHF_TreeCreator/tree_event_char", fname)
np_ev = tdf_ev.AsNumpy()
df_ev = pd.DataFrame(np_ev)
print (df_ev)