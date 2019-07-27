#!/usr/bin/env python

import pythia8
import pythiaext

def create_and_init_pythia(config_strings=[]):
	pythia = pythia8.Pythia()
	for s in config_strings:
		pythia.readString(s)
	for extra_s in ["Next:numberShowEvent = 0", "Next:numberShowInfo = 0", "Next:numberShowProcess = 0", "Next:numberCount = 0"]:
		pythia.readString(extra_s)
	if pythia.init():
		return pythia
	return None

def main():
	pythia = create_and_init_pythia(["PhaseSpace:pTHatMin = 2", "HardQCD:all = on"])
	sfoutname = "test_write_pythia_hepmc2.dat"
	pyhepmcwriter = pythiaext.Pythia8HepMCWrapper(sfoutname)
	for iEvent in range(100):
		if not pythia.next(): continue
		pyhepmcwriter.fillEvent(pythia)
	pythia.stat()
	print("[i] done writing to {}".format(sfoutname))

if __name__ == '__main__':
	main()
