#!/usr/bin/env python3

import os
import argparse
import ROOT as r
import math


def get_EPPS16_uncerts(flist, hname='hjetpt'):
	hset0 = None
	hdiffs = []
	for i, f in enumerate(flist):
		fin = r.TFile(f)
		hset = fin.Get(hname)
		if i == 0:
			hset0 = hset.Clone('hset0')
			hset0.SetDirectory(0)
			print('[i] set 0 set {}'.format(hset0))
			continue
		hsetdiff = hset.Clone('hset{}'.format(i))
		hsetdiff.SetDirectory(0)
		# print(hsetdiff)
		hsetdiff.Add(hset0, -1.)
		hdiffs.append(hsetdiff)
		fin.Close()

	hstat_plus = hset0.Clone('hstat_plus')
	hstat_minus = hset0.Clone('hstat_minus')

	for ib in range(1, hset0.GetNbinsX() + 1):
		diffplus = []
		diffminus = []
		for i,h in enumerate(hsetdiff):
			diff = hsetdiff.GetBinContent(ib)
			if i % 2 == 0:
				diffplus.append(diff)
			else:
				diffminus.append(diff)
		_dplus = sum([_d*_d for _d in diffplus])
		_dminus = sum([_d*_d for _d in diffminus])
		_dplus = math.sqrt(_dplus)
		_dminus = math.sqrt(_dminus)
		hstat_plus.SetBinContent(ib, hset0.GetBinContent(ib) + _dplus)
		hstat_minus.SetBinContent(ib, hset0.GetBinContent(ib) - _dminus)

	fout = r.TFile('epps16_uncerts.root', 'recreate')
	fout.cd()
	hstat_minus.Write()
	hstat_plus.Write()
	hset0.Write()
	fout.Close()


def main(args):
	fbasename = args.input
	flist = []
	for iset in range(41):
		fname = fbasename.replace('EPPS16_set0', 'EPPS16_set{0}'.format(iset))
		flist.append(fname)
	get_EPPS16_uncerts(flist, 'hjetpt')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jets with different with EPPS16nlo_CT14nlo_Pb208',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='root file with set0', type=str, required=True)
	args = parser.parse_args()
	main(args)