#!/usr/bin/env python

# copy that file somewhere and start generating...

from pythiautils import configuration as pyconf
import argparse
import os


def main(args):
	print(args)
	mycfg = []
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='some pythia simulation', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	args = parser.parse_args()	
	main(args)
