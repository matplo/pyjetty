""" eff_maker.py
Python script for converting the THnSparse of matches tracks to usable
histograms of efficiency and pT smearing (TODO) for use in fast simulation.
Author: Ezra Lesser (elesser@berkeley.edu)
"""

from __future__ import print_function

import ROOT


# Open ROOT files with relevant histograms
inf_true_name = "tracks_PhysPrim.root"
inf_det_name = "tracks.root"
inf_match_name = "tracks_Matched.root"

print("Opening ROOT files...")
inf_true = ROOT.TFile.Open(inf_true_name, "READ")
inf_match = ROOT.TFile.Open(inf_match_name, "READ")

print("Getting histograms...")
thn_true = inf_true.Get("tracks_PhysPrim;1")
thn_match = inf_match.Get("tracks_Matched;1")

# Tracking efficiency
print("Creating tracking efficiency plot...")
pt_true = thn_true.Projection(0, 3)
pt_true_findable = pt_true.ProjectionY("pt_true_findable", 2, 2)
pt_gen_matched = thn_match.Projection(0)
tr_eff = pt_gen_matched.Clone()
tr_eff.SetNameTitle("tr_eff", "tr_eff")
tr_eff.SetYTitle("Tracking Efficiency")
tr_eff.Divide(pt_gen_matched, pt_true_findable, 1, 1, "B")

# Save histograms to output file
print("Saving to output file...")
outf_name = "tr_eff.root"
outf = ROOT.TFile.Open(outf_name, "RECREATE")
outf.cd()
tr_eff.Write()
outf.Close()

# Close input files
inf_true.Close()
inf_match.Close()
