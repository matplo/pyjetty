#!/usr/bin/env python

from pythiautils import configuration as pyconf
import argparse
import os
import tqdm

def main(args):
	print(args)
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue
	pythia.stat()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='jet reco on alice data', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	
	main(args)
