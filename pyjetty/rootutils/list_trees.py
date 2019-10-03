#!/usr/bin/env python

import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os
import argparse
import pyjetty.rootutils as rootutils


def main():
	parser = argparse.ArgumentParser(description='list trees in a file', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--fname', help='path to a root file', default=None, type=str, required=True)
	parser.add_argument('-b', '--branches', help='list also branches', default=None, required=False, action="store_true")
	args = parser.parse_args()	
	with rootutils.QuietWarning():
		tu = rootutils.list_trees_dict(args.fname)
	for t in tu:
		print('Tree name: {}'.format(t.decode("utf-8")))
		if args.branches:
			_path = os.path.join(args.fname, t.decode("utf-8"))
			with rootutils.QuietWarning():
				tu = rootutils.TreeUtils(_path)
				bnames = tu.bnames
			for b in bnames:
				print ('  - {}'.format(b))


if __name__ == '__main__':
	main()