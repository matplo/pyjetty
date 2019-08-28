#!/usr/bin/env python3

import os
import argparse
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
import math
import dlist
import numpy
import copy

def make_ratio_graphs(fref, fepps16, hname):
	print ('[i] making [graphs] ratio of',fref, fepps16)
	f0 = r.TFile(fref)
	href = f0.Get(hname)
	grref = dlist.h_to_graph(href, False, True)
	f1 = r.TFile(fepps16)
	grepps0 = f1.Get("{}_set0_to_graph".format(hname))
	greppsU = f1.Get("{}_epps16_uncerts".format(hname))
	ryeh = []
	ryel = []
	exl = []
	exh = []
	ex = []
	ey = []
	for ib in range(greppsU.GetN()):
		ryeh.append(greppsU.GetErrorYhigh(ib))
		ryel.append(greppsU.GetErrorYlow(ib))
		exl.append(greppsU.GetErrorXlow(ib))
		exh.append(greppsU.GetErrorXhigh(ib))
		ex.append(grepps0.GetErrorX(ib))
		ey.append(grepps0.GetErrorY(ib))
	rx   = [_v for _v in greppsU.GetX()]
	ry   = [_v for _v in greppsU.GetY()]
	refy = [_v for _v in grref.GetY()]
	#print (ry, refy)
	for i in range(greppsU.GetN()):
		if refy[i] > 0:
			ry[i] 	= ry[i] / refy[i]
			ryel[i] = ryel[i] / refy[i]
			ryeh[i] = ryeh[i] / refy[i]
			ey[i]   = ey[i] / refy[i]
		else:
			ry[i] 	= 1.
			ryel[i] = 1.
			ryeh[i] = 1.
			ey[i]   = 1.

	foutname = fepps16.replace('uncerts_output', 'uncerts_ratio')
	fout = r.TFile(foutname, 'recreate')
	gr_ratio = dlist.make_graph_ae_xy('{}_epps16_uncerts_ratio'.format(hname), rx, ry, exl, exh, ryel, ryeh)
	gr_ratio.SetTitle('EPPS16 / default')
	gr_ratio.Write()
	gr_ratio0 = dlist.make_graph_xy('{}_epps16_uncerts_ratio_noErr'.format(hname), rx, ry, ex, ey)
	gr_ratio0.SetTitle('EPPS16 / default')
	gr_ratio0.Write()
	fout.Close()
	print ("[i] written", foutname)

def make_ratio(fref, fepps16):
	print ('[i] making ratio of',fref, fepps16)
	f0 = r.TFile(fref)
	f1 = r.TFile(fepps16)
	h0 = f0.Get('pt')
	hepps_set0 = f1.Get('set0')
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


def make_EPPS16_uncerts_new(flist, hname='hjetpt'):
	hset0 = None
	x = []
	y0 = []
	yA = []
	yeplus = []
	yeminus = []
	hset0 = None
	for i, f in enumerate(flist):
		if not os.path.isfile(f):
			print ('[w] missing', f)
			continue
		fin = r.TFile(f)
		hset = fin.Get(hname)
		if not hset:
			print ('[w] missing', hname, 'in', f)
			continue
		if i == 0:
			hset0 = hset.Clone('{}_set0'.format(hname))
			hset0.SetDirectory(0)
			print('[i] set 0 set {}'.format(hset0))
			for ibx in range(1, hset0.GetNbinsX()+1):
				if hset0.GetBinContent(ibx) > 0:
					x.append(hset0.GetBinCenter(ibx))
					y0.append(hset0.GetBinContent(ibx))
			yA.append(y0)
			continue
		if hset0 is None:
			print("[error] no set 0 found? - this will not work")
			return ""
		_ye = []
		for ibx in range(1, hset.GetNbinsX()+1):
			_ye.append(hset.GetBinContent(ibx))
		yA.append(_ye)

	# print(yA)
	nparams = int((len(flist) - 1) / 2)

	for ix in range(len(x)):
		dyplus = 0
		dyminus = 0
		for iset in range(0, nparams):
			_iset_plus = iset*2+1
			_iset_minus = iset*2+2
			# print(ix, _iset_plus, _iset_minus)
			_dyplus  = max(yA[_iset_plus][ix] - y0[ix], yA[_iset_minus][ix] - y0[ix], 0)
			_dyminus = max(y0[ix] - yA[_iset_plus][ix], y0[ix] - yA[_iset_minus][ix], 0)
			dyplus  += (_dyplus  * _dyplus)
			dyminus += (_dyminus * _dyminus)
		dyplus  = math.sqrt(dyplus)
		dyminus = math.sqrt(dyminus)
		yeplus.append(dyplus)
		yeminus.append(dyminus)

	foutname = os.path.join(flist[0].replace('EPPS16_set0' ,'uncerts'))
	fout = r.TFile(foutname, 'recreate')
	fout.cd()
	grset0 = dlist.h_to_graph(hset0, False, True)
	gr_ex = grset0.GetEX()
	gr_exl = [ex for ex in gr_ex]
	gr_exh = [ex for ex in gr_ex]
	grassymErr = dlist.make_graph_ae_xy('{}_epps16_uncerts'.format(hname), x, y0, gr_exl, gr_exh, yeminus, yeplus)
	hset0.Write()
	grset0.Write()
	grassymErr.Write()
	yplus   = [y0[i] + yeplus[i] for i in range(len(x))]
	grplus  = dlist.make_graph_xy('{}_epps16_uncerts_plus'.format(hname), x, yplus)
	grplus.Write()
	smallest_nonzero_float=numpy.nextafter(0,1)
	yminus  = [max(y0[i] - yeminus[i], 1e-20) for i in range(len(x))]
	grminus = dlist.make_graph_xy('{}_epps16_uncerts_minus'.format(hname), x, yminus)
	grminus.Write()
	fout.Close()
	print ("[i] written", foutname)
	return foutname


def make_EPPS16_uncerts(flist, hname='hjetpt'):
	hset0 = None
	hdiffs = []
	for i, f in enumerate(flist):
		if not os.path.isfile(f):
			print ('[w] missing', f)
			continue
		fin = r.TFile(f)
		hset = fin.Get(hname)
		if not hset:
			print ('[w] missing', hname, 'in', f)
			continue
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

	# print (len(hdiffs))
	for ib in range(1, hset0.GetNbinsX() + 1):
		diffplus = []
		diffminus = []
		for i,h in enumerate(hdiffs):
			diff = h.GetBinContent(ib)
			if i % 2 == 0:
				diffminus.append(diff)
			else:
				diffplus.append(diff)
			# print ('set', i)
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
	funcerts = make_EPPS16_uncerts_new(flist, 'pt')
	fref = os.path.join(fbasename.replace('EPPS16_set0_', ''))
	#make_ratio(fref, funcerts)
	make_ratio_graphs(fref, funcerts, 'pt')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jets with different with EPPS16nlo_CT14nlo_Pb208',
	                                 prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='root file with set0', type=str, required=True)
	args = parser.parse_args()
	main(args)
