#!/usr/bin/env python3

from __future__ import print_function
import tqdm
import argparse
import os
import numpy as np
import sys
import sys
import pyhepmc


def count_events(fname, hepmc_ver = 3):
	if hepmc_ver == 3:
		input_hepmc = pyhepmc.io.ReaderAscii(fname)
	if hepmc_ver == 2:
		input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(fname)

	if input_hepmc.failed():
		print ("[error] unable to read from {}".format(fname))
		sys.exit(1)

	pbar = tqdm.tqdm()
	event_hepmc = pyhepmc.GenEvent()
	while not input_hepmc.failed():
		ev = input_hepmc.read_event(event_hepmc)
		if input_hepmc.failed():
			break
		pbar.update(1)
  
	pbar.close()
	return pbar.n

def get_n_per_file(fname, hepmc_ver = 3):
	stat_file = fname 
	stat_file += '.ev_count'
	cached = False
	if os.path.exists(stat_file):
		with open(stat_file) as f:
			rlines = f.readlines()
		for l in rlines:
			if len(l) > 0:
				nev = int(l.split()[0])
				cached = True
				return nev, cached
	nev = count_events(fname, hepmc_ver)
	with open(stat_file, 'w') as f:
		f.writelines([ f'{nev}\n' ])
	return nev, cached


def main():
	parser = argparse.ArgumentParser(description='read hepmc and analyze eecs', prog=os.path.basename(__file__))
	parser.add_argument('-i', '--input', help='input file', default='low', type=str, required=True)
	parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
	args = parser.parse_args()	

	nev, cached = get_n_per_file(args.input, args.hepmc)
	print(f'[i] number of events (cached={cached}): {nev}')

if __name__ == '__main__':
	main()
