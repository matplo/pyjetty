#!/usr/bin/env python3

import os
import argparse
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
import math
import dlist

def make_ratio(fref, fepps16):
	print ('[i] making ratio of',fref, fepps16)
	f0 = r.TFile(fref)
	f1 = r.TFile(fepps16)
	h0 = f0.Get('hjetpt')
	hepps_set0 = f1.Get('hset0')
	hepps_plus = f1.Get('hstat_plus')
	hepps_minus = f1.Get('hstat_minus')

	h0.Sumw2()
	hepps_set0.Sumw2()
	hepps_plus.Sumw2()
	hepps_plus.Sumw2()

	dlist.scale_by_binwidth(h0)
	dlist.scale_by_binwidth(hepps_set0)
	dlist.scale_by_binwidth(hepps_plus)
	dlist.scale_by_binwidth(hepps_minus)

	hepps_set0.Divide(hepps_set0, h0, 1., 1., 'B')
	hepps_plus.Divide(hepps_plus, h0, 1., 1., 'B')
	hepps_minus.Divide(hepps_minus, h0, 1., 1., 'B')
	# hepps_set0.Divide(h0)
	# hepps_plus.Divide(h0)
	# hepps_minus.Divide(h0)

	hepps_set0.SetTitle('set0')
	hepps_plus.SetTitle('plus')
	hepps_minus.SetTitle('minus')

	foutname = fepps16.replace('uncerts_output', 'uncerts_ratio')
	fout = r.TFile(foutname, 'recreate')
	hepps_set0.Write()
	hepps_plus.Write()
	hepps_minus.Write()

	grset0 = dlist.h_to_graph(hepps_set0, False, True)
	gr_x = grset0.GetX()
	gr_y = grset0.GetY()
	gr_ex = grset0.GetEX()
	gr_exl = [ex for ex in gr_ex]
	gr_exh = [ex for ex in gr_ex]
	gr_eyl = [hepps_set0.GetBinContent(i) - hepps_minus.GetBinContent(i) for i in range(1, hepps_minus.GetNbinsX() + 1)]
	gr_eyh = [hepps_plus.GetBinContent(i) - hepps_set0.GetBinContent(i) for i in range(1, hepps_plus.GetNbinsX() + 1)]
	grassymErr = dlist.make_graph_ae_xy('epps16_uncerts_ratio', gr_x, gr_y, gr_exl, gr_exh, gr_eyl, gr_eyh)
	grassymErr.SetTitle('EPPS16 / default')
	grassymErr.Write()

	gr = dlist.make_graph_xy('epps16_uncerts_ratio_noErr', gr_x, gr_y)
	gr.SetTitle('EPPS16 / default')
	gr.Write()

	fout.Close()
	print ("[i] written", foutname)

def make_EPPS16_uncerts(flist, hname='hjetpt'):
	hset0 = None
	hdiffs = []
	for i, f in enumerate(flist):
		if not os.path.isfile(f):
			print ('[w] missing', f)
			continue
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

	grset0 = dlist.h_to_graph(hset0, False, True)

	# gr_x = [x*1.1 for x in grset0.GetX()]
	gr_x = grset0.GetX()
	gr_y = grset0.GetY()
	gr_ex = grset0.GetEX()
	gr_exl = [ex for ex in gr_ex]
	gr_exh = [ex for ex in gr_ex]
	gr_eyl = []
	gr_eyh = []

	print (len(hdiffs))
	for ib in range(1, hset0.GetNbinsX() + 1):
		diffplus = []
		diffminus = []
		for i,h in enumerate(hdiffs):
			diff = h.GetBinContent(ib)
			if i % 2 == 0:
				diffminus.append(diff)
			else:
				diffplus.append(diff)
			print ('set', i)
		_dplus = sum([_d*_d for _d in diffplus])
		_dminus = sum([_d*_d for _d in diffminus])
		_dplus = math.sqrt(_dplus)
		_dminus = math.sqrt(_dminus)
		print ('delta plus-minus', _dplus-_dminus)
		hstat_plus.SetBinContent(ib, hset0.GetBinContent(ib) + _dplus)
		hstat_minus.SetBinContent(ib, hset0.GetBinContent(ib) - _dminus)
		gr_eyl.append(_dminus)
		gr_eyh.append(_dplus)

	foutname = os.path.join(flist[0].replace('EPPS16_set0' ,'uncerts'))
	fout = r.TFile(foutname, 'recreate')
	fout.cd()
	hstat_minus.Write()
	hstat_plus.Write()
	hset0.Write()
	grset0.Write()
	grassymErr = dlist.make_graph_ae_xy('epps16_uncerts', gr_x, gr_y, gr_exl, gr_exh, gr_eyl, gr_eyh)
	grassymErr.Write()
	fout.Close()
	print ("[i] written", foutname)
	return foutname

def main(args):
	fbasename = args.input
	flist = []
	for iset in range(41):
		fname = fbasename.replace('EPPS16_set0', 'EPPS16_set{0}'.format(iset))
		flist.append(fname)
	funcerts = make_EPPS16_uncerts(flist, 'hjetpt')
	fref = os.path.join(fbasename.replace('EPPS16_set0_', ''))
	make_ratio(fref, funcerts)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jets with different with EPPS16nlo_CT14nlo_Pb208',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='root file with set0', type=str, required=True)
	args = parser.parse_args()
	main(args)
