#!/usr/bin/env python

import aleph
import aleph_utils

from tqdm import tqdm
import argparse
import os
import gzip

from pyjetty.mputils import pwarning, pinfo, perror, treewriter
from pyjetty.mputils import get_value

import ROOT
ROOT.gROOT.SetBatch()


def get_substr(l, t1, t2 = ' '):
	return l.split(t1)[1].split(t2)[0]


def get_event_info(l):
	run_number = 	get_value(get_substr(l, 'ALEPH_DATA RUN ='), 	int, -1)
	event_number = 	get_value(get_substr(l, 'EVENT '), 				int, -1)
	ecm = 			get_value(get_substr(l, 'ECM ='), 				int, -1)
	return run_number, event_number, ecm

def get_pvertex_info(l):
	vflag = 	get_value(get_substr(l, 'Primary vertex info flag ='), 	int, -1)
	vx = 	get_value(get_substr(l, 'vx ='), 	float, 0)
	vy = 	get_value(get_substr(l, 'vy ='), 	float, 0)
	ex = 	get_value(get_substr(l, 'ex ='), 	float, 0)
	ey = 	get_value(get_substr(l, 'ey =', '\n'), 	float, 0)
	return vflag, vx, vy, ex, ey

def get_part(l):
	px = get_value(get_substr(l, 'px='), 		float, 0)
	py = get_value(get_substr(l, 'py='), 		float, 0)
	pz = get_value(get_substr(l, 'pz='), 		float, 0)
	m  = get_value(get_substr(l, 'm='), 		float, 0)

	q  = get_value(get_substr(l, 'charge '), 	float, 0)
	pwflag  = get_value(get_substr(l, 'pwflag '), 	float, 0)
	d0  = get_value(get_substr(l, 'd0 '), 	float, 0)
	z0  = get_value(get_substr(l, 'z0 '), 	float, 0)
	ntpc  = get_value(get_substr(l, 'ntpc '), 	float, 0)
	nitc  = get_value(get_substr(l, 'nitc '), 	float, 0)
	nvdet  = get_value(get_substr(l, 'nvdet ', '\n'), 	float, 0)

	return px, py, pz, m, q, pwflag, d0, z0, ntpc, nitc, nvdet

def main(args):
	if args.output == 'default.root':
		args.output = args.input + '.root'

	pinfo('args', args)

	tw = treewriter.RTreeWriter(name = 'taleph', file_name = args.output)

	with gzip.open(args.input) as f:
		data = f.readlines()
		pinfo('number of lines read', len(data))
		for l in tqdm(data):
			lstr = l.decode("utf-8")
			if 'ALEPH_DATA RUN' in lstr:
				run_number, event_number, ecm = get_event_info(lstr)
				continue
			if 'Primary vertex info' in lstr:
				vflag, vx, vy, ex, ey = get_pvertex_info(lstr)
				continue
			if 'END_EVENT' in lstr:
				continue
			if 'px=' in lstr:
				px, py, pz, m, q, pwflag, d0, z0, ntpc, nitc, nvdet = get_part(lstr)
				tw.fill_branches(	run 	= run_number,
									event 	= event_number,
									ecm 	= ecm,
									vflag 	= vflag,
									vx 		= vx,
									vy 		= vy,
									ex 		= ex,
									ey 		= ey,
									px		= px,
									py		= py,
									pz		= pz,
									m		= m	,
									q		= q	,
									pwflag	= pwflag,
									d0		= d0,
									z0		= z0,
									ntpc	= ntpc,
									nitc	= nitc,
									nvdet	= nvdet)
				tw.fill_tree()
	tw.write_and_close()


def test_cxx(args):
	reader = aleph.Reader(args.input)
	nev = aleph_utils.get_n_events(args.input)
	for i in tqdm(range(nev)):
		if reader.read_next_event():
			e = reader.get_event()
		else:
			pinfo('no more events to read')
			break

def stream_particle(e, p, tw):
	# tlv = ROOT.TLorentzVector()
	# tlv.SetPxPyPzE(p.px(), p.py(), p.pz(), p.e())
	eh = e.get_header()
	tw.fill_branches(	
		run 	= eh.run(),
		event 	= eh.n(),
		ecm 	= eh.e(),
		vflag 	= eh.vflag(),
		vx 		= eh.vx(),
		vy 		= eh.vy(),
		ex 		= eh.ex(),
		ey 		= eh.ey(),
		px		= p.px(),
		py		= p.py(),
		pz		= p.pz(),
		m		= p.m(),
		e		= p.e(),
		pt		= p.pt(),
		# phi		= tlv.Phi(),
		# eta     = tlv.Eta(),
		q		= p.q(),
		pwflag	= p.pwflag(),
		d0		= p.d0(),
		z0		= p.z0(),
		ntpc	= p.ntpc(),
		nitc	= p.nitc(),
		nvdet	= p.nvdet())
	tw.fill_tree()


def test_cxx_gzip(args):

	pinfo(args)

	tw = treewriter.RTreeWriter(name = 'aleph', file_name = args.output)

	nev = aleph_utils.get_n_events_gzip(args.input)
	with gzip.open(args.input) as f:
		_data = f.readlines()
		data = aleph.StringVector()
		__ = [data.push_back("{}".format(s.decode("utf-8").strip('\n'))) for s in _data]
		pinfo(type(data))
		pinfo('number of lines read', len(data))
		reader = aleph.ReaderLines(data)
		for i in tqdm(range(nev)):
			if reader.read_next_event():
				e = reader.get_event()
				__ = [ stream_particle(e, p, tw) for p in e.get_particles()]
			else:
				pinfo('no more events to read')
				break
	tw.write_and_close()


def make_dict_row(e, p, dict_data):
	tlv = ROOT.TLorentzVector()
	tlv.SetPxPyPzE(p.px(), p.py(), p.pz(), p.e())
	eh = e.get_header()
	dict_data.append(
		{
		'run' : eh.run(),
		'event' : eh.n(),
		'ecm' : eh.e(),
		'vflag' : eh.vflag(),
		'vx' : eh.vx(),
		'vy' : eh.vy(),
		'ex' : eh.ex(),
		'ey' : eh.ey(),
		'px' : p.px(),
		'py' : p.py(),
		'pz' : p.pz(),
		'm' : p.m(),
		'e' : p.e(),
		'pt' : p.pt(),
		'phi' : tlv.Phi(),
		'eta' : tlv.Eta(),
		'q' : p.q(),
		'pwflag' : p.pwflag(),
		'd0' : p.d0(),
		'z0' : p.z0(),
		'ntpc' : p.ntpc(),
		'nitc' : p.nitc(),
		'nvdet' : p.nvdet()
		})


def write_csv(args):
	import csv
	csv_columns = [
		'run',
		'event',
		'ecm',
		'vflag',
		'vx',
		'vy',
		'ex',
		'ey',
		'px',
		'py',
		'pz',
		'm',
		'e',
		'pt',
		'phi',
		'eta',
		'q',
		'pwflag',
		'd0',
		'z0',
		'ntpc',
		'nitc',
		'nvdet']
	dict_data = []

	nev = aleph_utils.get_n_events_gzip(args.input)
	with gzip.open(args.input) as f:
		_data = f.readlines()
		data = aleph.StringVector()
		__ = [data.push_back("{}".format(s.decode("utf-8").strip('\n'))) for s in _data]
		pinfo(type(data))
		pinfo('number of lines read', data.size())
		reader = aleph.ReaderLines(data)
		for i in tqdm(range(nev)):
			if reader.read_next_event():
				e = reader.get_event()
				__ = [ make_dict_row(e, p, dict_data) for p in e.get_particles()]
			else:
				pinfo('no more events to read')
				break

	pinfo('writing csv file', args.output)
	csv_file = args.output
	try:
	    with open(csv_file, 'w') as csvfile:
	        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
	        writer.writeheader()
	        for data in dict_data:
	            writer.writerow(data)
	except IOError:
	    print("I/O error")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='ALEPH text to ROOT', prog=os.path.basename(__file__))
	required_args = parser.add_argument_group('required named arguments')
	required_args.add_argument('-i', '--input', help='input ALEPH gz file', required=True, type=str)
	parser.add_argument('-o', '--output', nargs='?', help='output file', default='default.root', type=str)
	parser.add_argument('--cxx', help='use cxx interface', action='store_true', default=False)
	args = parser.parse_args()
	if args.cxx:
		if '.gz' in args.input:
			if args.output[-4:] == '.csv':
				write_csv(args)
			else:
				test_cxx_gzip(args)
		else:
			test_cxx(args)
	else:
		main(args)