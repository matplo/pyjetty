#!/usr/bin/env python

import ROOT
import joblib
import pandas as pd
import mputils

def scale_by_bw(h):
	print(h.GetName())
	h.Sumw2()
	for ib in range(1, h.GetNbinsX()+1):
		v = h.GetBinContent(ib)
		h.SetBinContent(ib, v / h.GetBinWidth(ib))
		e = h.GetBinError(ib)
		h.SetBinError(ib, e / h.GetBinWidth(ib))

def plot_some(infile, outfile):
	data = joblib.load(infile)

	fout = ROOT.TFile(outfile, 'recreate')
	hpt = ROOT.TH1F('hpt', 'hpt', 10, mputils.logbins(1, 100, 11))
	[ hpt.Fill(x) for x in data['ptsub'] ]
	scale_by_bw(hpt)

	dsel = data.loc[(data['ptsub'] > 20) & (data['ptsub'] < 60)]

	hptsel = ROOT.TH1F('hptsel', 'hptsel', 10, mputils.logbins(1, 100, 11))
	[ hptsel.Fill(x) for x in dsel['ptsub'] ]
	scale_by_bw(hptsel)

	print(infile, len(dsel))

	#hsd0zg = ROOT.TH1F('hsd0zg', 'hsd0zg', 22, -1.5, 0.7)
	#hsd2zg = ROOT.TH1F('hsd2zg', 'hsd2zg', 22, -1.5, 0.7)

	hsd0zg = ROOT.TH1F('hsd0zg', 'hsd0zg', 12, 0, 0.6)
	[ hsd0zg.Fill(x) for x in dsel['sd0zg'] ]
	scale_by_bw(hsd0zg)
	hsd2zg = ROOT.TH1F('hsd2zg', 'hsd2zg', 12, 0, 0.6)
	[ hsd2zg.Fill(x) for x in dsel['sd2zg'] ]
	scale_by_bw(hsd2zg)

	hsd0Rg = ROOT.TH1F('hsd0Rg', 'hsd0Rg', 20, 0, 0.6)
	[ hsd0Rg.Fill(x) for x in dsel['sd0Rg'] ]
	scale_by_bw(hsd0Rg)
	hsd2Rg = ROOT.TH1F('hsd2Rg', 'hsd2Rg', 20, 0, 0.6)
	[ hsd2Rg.Fill(x) for x in dsel['sd2Rg'] ]
	scale_by_bw(hsd2Rg)

	hsd0thg = ROOT.TH1F('hsd0thg', 'hsd0thg', 15, 0, 1.5)
	[ hsd0thg.Fill(x) for x in dsel['sd0thetag'] ]
	scale_by_bw(hsd0thg)
	hsd2thg = ROOT.TH1F('hsd2thg', 'hsd2thg', 15, 0, 1.5)
	[ hsd2thg.Fill(x) for x in dsel['sd2thetag'] ]
	scale_by_bw(hsd2thg)

	fout.Write()
	fout.Close()

def main():
	plot_some('merged_rg_PbPb_std.pd', "hist_std.root")
	plot_some('merged_rg_PbPb_cs004.pd', "hist_cs004.root")
	plot_some('merged_rg_PbPb_cs404.pd', "hist_cs404.root")

if __name__ == '__main__':
	main()
