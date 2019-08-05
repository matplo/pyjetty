#!/usr/bin/env python

from __future__ import print_function

import tqdm
import argparse
import os
import numpy as np

import pythia8

from pythiautils import configuration as pyconf

parser = argparse.ArgumentParser(description='jet reco on alice data', prog=os.path.basename(__file__))
pyconf.add_standard_pythia_args(parser)
args = parser.parse_args()	

mycfg = ['PhaseSpace:pThatMin = 10']
pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)
if args.nev < 10:
	args.nev = 10
for i in tqdm.tqdm(range(args.nev)):
	if not pythia.next():
		continue

pythia.stat()
